{- Calculate contig and scaffold contiguity and correctness metrics.
 - Written by Shaun Jackman.
 -}

import Control.Arrow ((***))
import Control.Monad (when)
import Data.Bits ((.&.))
import Data.ByteString.Lazy.Char8 (ByteString)
import qualified Data.ByteString.Lazy.Char8 as S
import Data.Char (isDigit)
import Data.Function (on)
import Data.List (find, findIndex, groupBy, intercalate, partition, sortBy, span)
import Data.Maybe (fromJust, fromMaybe)
import System.Console.GetOpt
import System.Environment (getArgs)
import System.Exit (exitFailure, exitSuccess)

-- Return groups of equivalent elements. Unlike Data.List.groupBy,
-- consecutive pairs of elements are compared, and the relation need
-- not be transitive.
groupBy' rel [] = []
groupBy' rel (x:xs) = (x:ys) : groupBy' rel zs
  where
  (ys,zs) = groupByAux x xs
  groupByAux x0 (x:xs) | rel x0 x = (x:ys, zs)
    where (ys,zs) = groupByAux x xs
  groupByAux y xs = ([], xs)

-- Return pairs of elements from a list.
pairs :: [a] -> [(a, a)]
pairs [] = []
pairs (x:y:ys) = (x, y) : pairs ys

-- Parse an integer from a ByteString.
readS :: ByteString -> Int
readS s = x where Just (x, _) = S.readInt s

-- Return the rank of the smallest element x for which the sum of the elements x or larger is at least c.
rankSumAtLeast :: Int -> [Int] -> Int
rankSumAtLeast c xs = i
  where
  i = fromMaybe 0 $ findIndex ((>= c) . snd) $ zip ws $ scanl1 (+) ws
  ws = sortBy (flip compare) xs

-- Return the smallest element x for which the sum of the elements x or larger is at least c.
sumAtLeast :: Int -> [Int] -> Int
sumAtLeast c xs = x
  where
  (x, _) = fromMaybe (0, 0) $ find ((>= c) . snd) $ zip ws $ scanl1 (+) ws
  ws = sortBy (flip compare) xs

-- Calculate NGx.
ngx :: Double -> Int -> [Int] -> Int
ngx x g = sumAtLeast $ ceiling $ x * fromIntegral g

-- Calculate N50.
n50 :: [Int] -> Int
n50 xs = ngx 0.5 (sum xs) xs

-- Calculate LGx.
lgx :: Double -> Int -> [Int] -> Int
lgx x g = rankSumAtLeast $ ceiling $ x * fromIntegral g

-- Calculate L50.
l50 :: [Int] -> Int
l50 xs = lgx 0.5 (sum xs) xs

-- A SAM record.
-- qname flag rname pos mapq cigar rnext pnext tlen seq qual
data SAM = SAM {
  qname :: ByteString,
  flag :: Int,
  rname :: ByteString,
  pos :: Int,
  mapq :: Int,
  cigar :: ByteString,
  --rnext :: ByteString,
  --pnext :: Int,
  --tlen :: Int,
  --qseq :: ByteString,
  --qual :: ByteString
  seqLength :: Int
}

-- Return whether the unmapped flag is set.
isUnmapped :: SAM -> Bool
isUnmapped x = flag x .&. 4 /= 0

-- Return whether the reverse complement flag is set.
isRC :: SAM -> Bool
isRC x = flag x .&. 16 /= 0

-- Parse the CIGAR string.
readCigar :: ByteString -> [(Int, Char)]
readCigar s = map (readS *** S.head) . pairs
  . S.groupBy ((==) `on` isDigit) $ s

-- Return the left and right soft-clipping.
getSoftClip :: SAM -> (Int, Int)
getSoftClip sam = (
    if elem (snd x) "HS" then fst x else 0,
    if elem (snd y) "HS" then fst y else 0)
  where
  (x, y) = (head xs, last xs)
  xs = readCigar $ cigar sam

-- Return the sum of the specified CIGAR elements.
cigarLength :: [Char] -> ByteString -> Int
cigarLength ops s = sum [n | (n, op) <- readCigar s, op `elem` ops]

-- Return the length of the query alignment.
qLength :: SAM -> Int
qLength = cigarLength "IM" . cigar

-- Return the start position of the query alignment oriented to agree
-- with the target.
qStart' = fst . getSoftClip

-- Return the end position of the query alignment oriented to agree
-- with the target.
qEnd' x = qStart' x + qLength x

-- Return the start position of the query alignment.
qStart :: SAM -> Int
qStart x = (if isRC x then snd else fst) $ getSoftClip x

-- Return the end position of the query alignment.
qEnd :: SAM -> Int
qEnd x = qStart x + qLength x

-- Return the length of the target alignment.
tLength :: SAM -> Int
tLength = cigarLength "DMN" . cigar

-- Return the start of the target alignment.
tStart :: SAM -> Int
tStart = pos

-- Return the end position of the target alignment.
tEnd :: SAM -> Int
tEnd x = tStart x + tLength x

-- Return the start position of the target alignment oriented to agree
-- with the query. The coordinate system is [-length, 0).
tStart' x = if isRC x then 0 - tEnd x else tStart x

-- Return the end position of the target alignment oriented to agree
-- with the query. The coordinate system is [-length, 0).
tEnd' x = tStart' x + tLength x

-- Parse a SAM record.
readSAM :: ByteString -> SAM
readSAM s = SAM qname (readS flag) rname (readS pos)
    (readS mapq) cigar lengthQseq
  where
  lengthQseq = if S.head cigar == '*'
    then fromIntegral $ S.length qseq
    else cigarLength "HIMS" cigar
  (qname:flag:rname:pos:mapq:cigar:_:_:_:qseq:_) = S.words s

-- Print a SAM record.
showSAM :: SAM -> String
showSAM (SAM qname flag rname pos mapq cigar _) = intercalate "\t"
    [S.unpack qname, show flag, S.unpack rname, show pos,
      show mapq, S.unpack cigar, "*", "0", "0", "*", "*"]

-- Exclude alignments that overlap a long alignment by 50%.
excludeOverlaps :: [SAM] -> [SAM]
excludeOverlaps xs = reverse . foldl accum [] $ xs
  where
  accum ys x = if any overlapx ys then ys else x:ys
      where
      overlapx y = end - start > (qLength x) `div` 2
          where
          start = max (qStart x) (qStart y)
          end = min (qEnd x) (qEnd y)

-- Compare the target position.
compareTStart :: SAM -> SAM -> Ordering
compareTStart a b = compare (rname a, pos a) (rname b, pos b)

-- Patch gaps in the alignments that are shorter than 500 bp.
patchGaps :: Int -> [SAM] -> [SAM]
patchGaps optPatchGaps ws = reverse . foldl accum [] $ xs
  where
  xs = sortBy compareTStart ws
  accum [] x = [x]
  accum (y:ys) x = if isSmallGap y x then x:ys else x:y:ys
  isSmallGap q p = (rname p, isRC p) == (rname q, isRC q)
      && max gapt gapq < optPatchGaps
    where
    gapt = abs (tStart p - tEnd q)
    gapq = abs (qStart' p - qEnd' q)

-- Return whether the pair of aligned contigs are colinear.
isColinear :: SAM -> SAM -> Bool
isColinear a b = (rname a, isRC a) == (rname b, isRC b)
  && tStart' a < tStart' b && tEnd' a < tEnd' b

-- Return pairs of non-colinear alignments.
filterNonColinear :: [SAM] -> [(SAM, SAM)]
filterNonColinear xs = filter (not . uncurry isColinear)
  $ zip xs (tail xs)

-- The command line options.
data Opt = OptAlignmentLength Int | OptContigLength Int | OptGenomeSize Int | OptMapq Int | OptPatchGaps Int |
  OptSAM | OptSAMContigs | OptSAMScaffolds | OptText | OptTSV |
  OptHelp | OptVersion
  deriving Eq
options :: [OptDescr Opt]
options = [
  Option ['a'] ["alignment-length"] (ReqArg (OptAlignmentLength . read) "N")
    "exclude alignments with aligned query length shorter than N bp [500]",
  Option ['l'] ["contig-length"] (ReqArg (OptContigLength . read) "N")
    "exclude contigs shorter than N bp [200]",
  Option ['G'] ["genome-size"] (ReqArg (OptGenomeSize . read) "N")
    "expected genome size used to calculate NG50 [0]",
  Option ['q'] ["mapq"] (ReqArg (OptMapq . read) "N")
    "exclude alignments with mapq less than N [10]",
  Option ['g'] ["patch-gaps"] (ReqArg (OptPatchGaps . read) "N")
    "join alignments separated by a gap shorter than INT bp (with --text only) [500]",
  Option ['p'] ["sam"] (NoArg OptSAM)
    "output contig and scaffold breakpoints in SAM format",
  Option [] ["sam-contigs"] (NoArg OptSAMContigs)
    "output contig breakpoints in SAM format",
  Option [] ["sam-scaffolds"] (NoArg OptSAMScaffolds)
    "output scaffold breakpoints in SAM format",
  Option [] ["text"] (NoArg OptText)
    "output report in plain text format",
  Option [] ["tsv"] (NoArg OptTSV)
    "output report in TSV format [default]",
  Option [] ["help"] (NoArg OptHelp)
    "display this help and exit",
  Option [] ["version"] (NoArg OptVersion)
    "display version information and exit" ]
data Format = FormatSAM | FormatSAMContigs | FormatSAMScaffolds | FormatText | FormatTSV deriving Eq
data Options = Options {
  optAlignmentLength :: Int,
  optContigLength :: Int,
  optGenomeSize :: Int,
  optMapq :: Int,
  optPatchGaps :: Int,
  optFormat :: Format
}
defaultOptions = Options {
  optAlignmentLength = 500,
  optContigLength = 200,
  optGenomeSize = 0,
  optMapq = 10,
  optPatchGaps = 500,
  optFormat = FormatTSV
}

-- Parse the command line options.
parseOptions :: [Opt] -> Options
parseOptions = foldl parseOption defaultOptions
  where
  parseOption opt x = case x of
    OptAlignmentLength n -> opt { optAlignmentLength = n }
    OptContigLength n -> opt { optContigLength = n }
    OptGenomeSize n -> opt { optGenomeSize = n }
    OptMapq n -> opt { optMapq = n }
    OptPatchGaps n -> opt { optPatchGaps = n }
    OptSAM -> opt { optFormat = FormatSAM }
    OptSAMContigs -> opt { optFormat = FormatSAMContigs }
    OptSAMScaffolds -> opt { optFormat = FormatSAMScaffolds }
    OptText -> opt { optFormat = FormatText }
    OptTSV -> opt { optFormat = FormatTSV }

-- Parse the command line arguments.
parseArgs :: IO (Options, [String])
parseArgs = do
  args <- getArgs
  case getOpt Permute options args of
    (opts, files, []) -> if OptHelp `elem` opts then help else
      if OptVersion `elem` opts then putStr version >> exitSuccess else
      return (parseOptions opts, files)
    (_, _, errs) -> error (concat errs ++ tryHelp)
  where
  help = putStr (usageInfo usage options) >> exitSuccess
  tryHelp = "Try 'abyss-samtobreak --help' for more information."
  version = "abyss-samtobreak (ABySS) 2.1.3\n"
  usage = "Usage: samtobreak [OPTION]... [FILE]...\n\
\Calculate contig and scaffold contiguity and correctness metrics.\n"

-- Calculate contig and scaffold contiguity and correctness metrics.
printStats :: Integer -> FilePath -> Options -> IO ()
printStats recordIndex path (Options optAlignmentLength optContigLength optGenomeSize optMapq optPatchGaps optFormat) = do
  s <- S.readFile path
  when (S.null s) $ error $ "`" ++ path ++ "' is empty"

  let
    -- Parse the SAM file and discard short contigs.
    isHeader x = S.head x == '@'
    (headers, alignments) = span isHeader . S.lines $ s
    allContigs = map readSAM $ alignments

    -- Calculate the size of the reference genome.
    sq = filter ((S.pack "@SQ" ==) . head) . map S.words $ headers
    referenceLengths = map (readS . S.drop 3 . fromJust . find (S.isPrefixOf $ S.pack "LN:")) sq
    reference_bases = sum referenceLengths

    -- Keep the first (primary) alignment of each contig.
    primaryContigs = map head . groupBy ((==) `on` qname) $ allContigs

    -- Group the contigs into scaffolds by their name.
    isScaffoldDelimiter c = c == '_' || c == ':'
    scaffoldName = S.takeWhile (not . isScaffoldDelimiter) . qname
    primaryScaffolds = groupBy ((==) `on` scaffoldName) $ primaryContigs

    -- Separate the mapped and unmapped contigs, and discard short contigs.
    isLong x = seqLength x >= optContigLength
    (unmapped, mapped) = partition isUnmapped $ filter isLong allContigs

    -- Exclude overlapping alignments.
    excluded = map excludeOverlaps
      . groupBy ((==) `on` qname) $ mapped
    concatExcluded = concat excluded
    alignedContigLengths = map qLength concatExcluded

    -- Keep long alignments with high mapping quality.
    isGood x = mapq x >= optMapq && qLength x >= optAlignmentLength
    good = filter (not . null) . map (filter isGood) $ excluded

    -- Discard contigs that have more than one alignment.
    -- Group contigs into scaffolds by their name.
    oneHit = concat . filter ((== 1) . length) $ good
    scaffs = groupBy ((==) `on` scaffoldName) $ oneHit

  -- Print SAM headers.
  when (optFormat == FormatSAM || optFormat == FormatSAMContigs || optFormat == FormatSAMScaffolds) (do
    S.putStr $ S.unlines headers)

  -- Print contig breakpoints.
  when (optFormat == FormatSAMContigs || optFormat == FormatSAM) (do
    mapM_ (putStrLn . showSAM) $ concat $ filter ((> 1) . length) $ good
    )

  -- Print scaffold breakpoints.
  when (optFormat == FormatSAMScaffolds || optFormat == FormatSAM) (do
    putStr $ concat
      $ map (\(a, b) -> unlines [showSAM a, showSAM b])
      $ concatMap filterNonColinear scaffs
    )

  -- Exit after printing the breakpoints.
  when (optFormat == FormatSAM || optFormat == FormatSAMContigs || optFormat == FormatSAMScaffolds) (do
    exitSuccess)

  let
    genomeSize = if optGenomeSize > 0 then optGenomeSize else reference_bases
    ng50 = ngx 0.5 genomeSize
    lg50 = lgx 0.5 genomeSize

    -- Contig metrics
    mapped_contigs = length good
    unmapped_contigs = length unmapped
    mapped_bases = sum alignedContigLengths
    unmapped_contig_bases = sum . map seqLength $ unmapped
    contigLengths = map seqLength $ filter isLong primaryContigs
    contig_n50 = n50 contigLengths
    contig_na50 = n50 alignedContigLengths
    contig_ng50 = ng50 contigLengths
    contig_nga50 = ng50 alignedContigLengths
    contig_max = maximum contigLengths
    contig_aligned_max = maximum alignedContigLengths
    contig_l50 = l50 contigLengths
    contig_la50 = l50 alignedContigLengths
    contig_lg50 = lg50 contigLengths
    contig_lga50 = lg50 alignedContigLengths
    contigs = length contigLengths
    contig_breakpoints = length (concat good) - length good

    -- Scaffold metrics
    scaffoldLengths = filter (>= optContigLength) . map (sum . map seqLength) $ primaryScaffolds
    colinearScaffs = concatMap (groupBy' isColinear) scaffs
    alignedScaffoldLengths = map (sum . map qLength) $ colinearScaffs
    scaffold_n50 = n50 scaffoldLengths
    scaffold_na50 = n50 alignedScaffoldLengths
    scaffold_ng50 = ng50 scaffoldLengths
    scaffold_nga50 = ng50 alignedScaffoldLengths
    scaffold_max = maximum scaffoldLengths
    scaffold_aligned_max = maximum alignedScaffoldLengths
    scaffold_l50 = l50 scaffoldLengths
    scaffold_la50 = l50 alignedScaffoldLengths
    scaffold_lg50 = lg50 scaffoldLengths
    scaffold_lga50 = lg50 alignedScaffoldLengths
    scaffolds = length scaffoldLengths
    total_bases = sum scaffoldLengths
    scaffold_breakpoints = length colinearScaffs - length scaffs
    total_breakpoints = contig_breakpoints + scaffold_breakpoints

  when (optFormat == FormatText) (do
    when (recordIndex > 0) (putChar '\n')

    putStrLn $ "File: " ++ path

    putStr "Number of unmapped contigs: "
    print unmapped_contigs

    putStr "Total length of unmapped contigs: "
    print unmapped_contig_bases

    putStr "Number of alignments dropped due to excessive overlaps: "
    print $ length mapped - length concatExcluded

    putStr "Mapped contig bases: "
    print mapped_bases

    putStr $ "Contig N50: "
    print contig_n50

    putStr $ "Contig NA50: "
    print contig_na50

    putStr $ "Contig NG50: "
    print contig_ng50

    putStr $ "Contig NGA50: "
    print contig_nga50

    putStr $ "Contig max: "
    print contig_max

    putStr $ "Contig aligned max: "
    print contig_aligned_max

    putStr "Number of break points: "
    print $ length concatExcluded - length excluded

    putStr $ "Number of Q10 break points longer than " ++ show optAlignmentLength ++ " bp: "
    print contig_breakpoints

    -- Patch small gaps.
    when (optPatchGaps > 0) (do
      let
        patched = map (patchGaps optPatchGaps) excluded
        concatPatched = concat patched

        patchedGood = filter (not . null) . map (filter isGood) $ patched
        concatPatchedGood = concat patchedGood

      putStr "Number of break points after patching gaps shorter than 500 bp: "
      print $ length concatPatched - length patched

      putStr "Number of Q10 break points longer than 500 bp after gap patching: "
      print $ length concatPatchedGood - length patchedGood
      )

    putStr $ "Scaffold N50: "
    print scaffold_n50

    putStr $ "Scaffold NA50: "
    print scaffold_na50

    putStr $ "Scaffold NG50: "
    print scaffold_ng50

    putStr $ "Scaffold NGA50: "
    print scaffold_nga50

    putStr $ "Scaffold max: "
    print scaffold_max

    putStr $ "Scaffold aligned max: "
    print scaffold_aligned_max

    putStr $ "Number of Q10 scaffold breakpoints longer than " ++ show optAlignmentLength ++ " bp: "
    print scaffold_breakpoints

    putStr "Number of contig and scaffold breakpoints: "
    print total_breakpoints
    )

  when (optFormat == FormatTSV) (do
    when (recordIndex == 0) (
      putStrLn ("File"
        ++ "\tContig_N50"
        ++ "\tContig_NA50"
        ++ "\tContig_NG50"
        ++ "\tContig_NGA50"
        ++ "\tContig_max"
        ++ "\tContig_aligned_max"
        ++ "\tScaffold_N50"
        ++ "\tScaffold_NA50"
        ++ "\tScaffold_NG50"
        ++ "\tScaffold_NGA50"
        ++ "\tScaffold_max"
        ++ "\tScaffold_aligned_max"
        ++ "\tContig_breakpoints"
        ++ "\tScaffold_breakpoints"
        ++ "\tTotal_breakpoints"
        ++ "\tContig_L50"
        ++ "\tContig_LA50"
        ++ "\tContig_LG50"
        ++ "\tContig_LGA50"
        ++ "\tScaffold_L50"
        ++ "\tScaffold_LA50"
        ++ "\tScaffold_LG50"
        ++ "\tScaffold_LGA50"
        ++ "\tReference_bases"
        ++ "\tTotal_bases"
        ++ "\tMapped_bases"
        ++ "\tUnmapped_contig_bases"
        ++ "\tContigs"
        ++ "\tMapped_contigs"
        ++ "\tUnmapped_contigs"
        ++ "\tScaffolds"
        ++ "\talignment-length"
        ++ "\tcontig-length"
        ++ "\tmapq"
        ++ "\tgenome-size")
      )
    putStr path
    putStr $ '\t' : show contig_n50
    putStr $ '\t' : show contig_na50
    putStr $ '\t' : show contig_ng50
    putStr $ '\t' : show contig_nga50
    putStr $ '\t' : show contig_max
    putStr $ '\t' : show contig_aligned_max
    putStr $ '\t' : show scaffold_n50
    putStr $ '\t' : show scaffold_na50
    putStr $ '\t' : show scaffold_ng50
    putStr $ '\t' : show scaffold_nga50
    putStr $ '\t' : show scaffold_max
    putStr $ '\t' : show scaffold_aligned_max
    putStr $ '\t' : show contig_breakpoints
    putStr $ '\t' : show scaffold_breakpoints
    putStr $ '\t' : show total_breakpoints
    putStr $ '\t' : show contig_l50
    putStr $ '\t' : show contig_la50
    putStr $ '\t' : show contig_lg50
    putStr $ '\t' : show contig_lga50
    putStr $ '\t' : show scaffold_l50
    putStr $ '\t' : show scaffold_la50
    putStr $ '\t' : show scaffold_lg50
    putStr $ '\t' : show scaffold_lga50
    putStr $ '\t' : show reference_bases
    putStr $ '\t' : show total_bases
    putStr $ '\t' : show mapped_bases
    putStr $ '\t' : show unmapped_contig_bases
    putStr $ '\t' : show contigs
    putStr $ '\t' : show mapped_contigs
    putStr $ '\t' : show unmapped_contigs
    putStr $ '\t' : show scaffolds
    putStr $ '\t' : show optAlignmentLength
    putStr $ '\t' : show optContigLength
    putStr $ '\t' : show optMapq
    putStr $ '\t' : show genomeSize
    putChar '\n'
    )

-- Calculate contig and scaffold contiguity and correctness metrics.
main :: IO ()
main = do
  (opt, files) <- parseArgs
  case files of
    [] -> printStats 0 "/dev/stdin" opt
    otherwise -> mapM_ (\(recordIndex, path) -> do
      printStats recordIndex path opt) (zip [0..] files)
