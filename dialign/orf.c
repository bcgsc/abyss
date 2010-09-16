/**
 *
 * orf.c: find the longest Open Reading Frame  
 *
 * 2004-08-24  Dorothea Emig Volker Menrad
 *             
 */

#include <stdlib.h>
#include <stdio.h>

#include "parameters.h"
#include "struct.h"
#include "translate.h"
#include "orf.h"
#include "io.h"

struct seq_col* set_longest_orf(struct seq_col *in_seq_col)
{
	int i;
	
	struct seq_col *ret_seq_col;
	struct seq *seq, *input_seq;
	if((ret_seq_col = calloc(1,sizeof(struct seq_col)))==NULL)
	{
		error("set_longest_orf(): Out of memory ");
	}
	if((ret_seq_col->seqs = calloc((in_seq_col->length),sizeof(struct seq)))==NULL)
	{
		error("set_longest_orf(): Out of memory ");
	}

	input_seq = &(in_seq_col->seqs[0]);
	for(i=0; i!= in_seq_col->length; ++i)
	{
		seq = orf_finder(input_seq);
		++input_seq;
		ret_seq_col->seqs[i].dna_num = seq->dna_num;
		ret_seq_col->seqs[i].orf_frame = seq->orf_frame;
		ret_seq_col->seqs[i].data = seq->data;
		ret_seq_col->seqs[i].name = seq->name;
		ret_seq_col->seqs[i].num = seq->num;
		ret_seq_col->seqs[i].length = seq->length;
		++seq;
		
	}
	ret_seq_col->length=in_seq_col->length;
	
	return ret_seq_col;
}


struct seq* orf_finder(struct seq *in_seq)
{
	char flag_no_start = 1; 
	struct seq *ret_seq;
	struct orf **orfs;
	int i, max_length, max_orf;
	char a,b,c,d,e,f; // first three possible reading frames and last three
	int help, j, x;                                           //          XXX........XXX
	j = (in_seq->length)/3;
	help = in_seq->length-1; // index of last element in sequence
	a=0;b=0;c=0;d=0;e=0;f=0;


	if((orfs = (calloc (12, sizeof(struct orf *)))) == NULL) // We need 12 because we look for the longest orf in every reading frame
															 // and save the longest for every reading frame = 6 
															 // orf[0] && orf[1] for reading frame a
															// orf[2] && orf[3] for reading frame b	
														    // orf[4] && orf[5] for reading frame c
															// orf[6] && orf[7] for reading frame d	
															// orf[8] && orf[9] for reading frame e
															// orf[10] && orf[11] for reading frame f	
	{														

		error("orf_finder(): Out of memory !");
	}

	for (x = 0; x != 12; ++x)
	{
		if((orfs[x] = (calloc(1,sizeof(struct orf))))==NULL)
		{
		error("orf_finder(): Out of memory !");
		}
		if((orfs[x]->sequence = (calloc(j+1,sizeof(char))))==NULL)
		{
		error("orf_finder(2): Out of memory !");
		}
		if((orfs[x]->dna_num = (calloc(j+1,sizeof(char))))==NULL)
		{
		error("orf_finder(3): Out of memory !");
		}

		orfs[x]->finish = 0;
	}

	for (i = 0; i < (j*3)-3; ++i)
	{
/****************************************************************************************/
		if (a == 1)
		{
			if(!(orfs[0]->finish))
			{
				if((in_seq->data[i] == 'T' && in_seq->data[i+1] == 'A' && in_seq->data[i+2] == 'A') ||
			 (in_seq->data[i] == 'U' && in_seq->data[i+1] == 'A' && in_seq->data[i+2] == 'A') ||	 
			(in_seq->data[i] == 'T' && in_seq->data[i+1] == 'A' && in_seq->data[i+2] == 'G') ||
			 (in_seq->data[i] == 'U' && in_seq->data[i+1] == 'A' && in_seq->data[i+2] == 'G') ||
			(in_seq->data[i] == 'T' && in_seq->data[i+1] == 'G' && in_seq->data[i+2] == 'A') ||
			 (in_seq->data[i] == 'U' && in_seq->data[i+1] == 'G' && in_seq->data[i+2] == 'A') ) //Stop-Codon then end orf
				{
					orfs[0]->finish = 1; // end orf
					a=0; // wait for new ATG-Start-Codon
				}
					orfs[0]->sequence[(orfs[0]->length)] = translate(in_seq->data[i], in_seq->data[i+1], in_seq->data[i+2], &(orfs[0]->dna_num[(orfs[0]->length)]) );	// translate in aminoacid and safe in orf
					++orfs[0]->length; 
			}
			else 
			{
				if((in_seq->data[i] == 'T' && in_seq->data[i+1] == 'A' && in_seq->data[i+2] == 'A') ||
			 (in_seq->data[i] == 'U' && in_seq->data[i+1] == 'A' && in_seq->data[i+2] == 'A') ||
			(in_seq->data[i] == 'T' && in_seq->data[i+1] == 'A' && in_seq->data[i+2] == 'G') ||
			 (in_seq->data[i] == 'U' && in_seq->data[i+1] == 'A' && in_seq->data[i+2] == 'G') ||
			(in_seq->data[i] == 'T' && in_seq->data[i+1] == 'G' && in_seq->data[i+2] == 'A') ||
			 (in_seq->data[i] == 'U' && in_seq->data[i+1] == 'G' && in_seq->data[i+2] == 'A') )
				{
					orfs[1]->finish = 1;
					a=0;
				}
				orfs[1]->sequence[(orfs[1]->length)] = translate(in_seq->data[i], in_seq->data[i+1], in_seq->data[i+2],&(orfs[1]->dna_num[(orfs[1]->length)]));	
				++orfs[1]->length;
			}		
		}
		
/****************************************************************************************/	

		if (b == 1)
		{
			if(!(orfs[2]->finish))
			{
				if((in_seq->data[i+1] == 'T' && in_seq->data[i+2] == 'A' && in_seq->data[i+3] == 'A') ||
			 (in_seq->data[i+1] == 'U' && in_seq->data[i+2] == 'A' && in_seq->data[i+3] == 'A') ||	 
			(in_seq->data[i+1] == 'T' && in_seq->data[i+2] == 'A' && in_seq->data[i+3] == 'G') ||
			 (in_seq->data[i+1] == 'U' && in_seq->data[i+2] == 'A' && in_seq->data[i+3] == 'G') ||
			(in_seq->data[i+1] == 'T' && in_seq->data[i+2] == 'G' && in_seq->data[i+3] == 'A') ||
			 (in_seq->data[i+1] == 'U' && in_seq->data[i+2] == 'G' && in_seq->data[i+3] == 'A') ) //Stop-Codon then end orf
				{
					orfs[2]->finish = 1; // end orf
					b=0; // wait for new ATG-Start-Codon
				}
				orfs[2]->sequence[(orfs[2]->length)] = translate(in_seq->data[i+1], in_seq->data[i+2], in_seq->data[i+3],&(orfs[2]->dna_num[(orfs[2]->length)]) );	// translate in aminoacid and safe in orf
				++orfs[2]->length; 
			}
			else 
			{
				if((in_seq->data[i+1] == 'T' && in_seq->data[i+2] == 'A' && in_seq->data[i+3] == 'A') ||
			 (in_seq->data[i+1] == 'U' && in_seq->data[i+2] == 'A' && in_seq->data[i+3] == 'A') ||
			(in_seq->data[i+1] == 'T' && in_seq->data[i+2] == 'A' && in_seq->data[i+3] == 'G') ||
			 (in_seq->data[i+1] == 'U' && in_seq->data[i+2] == 'A' && in_seq->data[i+3] == 'G') ||
			(in_seq->data[i+1] == 'T' && in_seq->data[i+2] == 'G' && in_seq->data[i+3] == 'A') ||
			 (in_seq->data[i+1] == 'U' && in_seq->data[i+2] == 'G' && in_seq->data[i+3] == 'A') )
				{
					orfs[3]->finish = 1;
					b=0;
				}
				orfs[3]->sequence[(orfs[3]->length)] = translate(in_seq->data[i+1], in_seq->data[i+2], in_seq->data[i+3],&(orfs[3]->dna_num[(orfs[3]->length)]) );	
				++orfs[3]->length;
			}		
		}

/****************************************************************************************/	

		if (c == 1)
		{
			if(!(orfs[4]->finish))

			{
				if((in_seq->data[i+2] == 'T' && in_seq->data[i+3] == 'A' && in_seq->data[i+4] == 'A') ||
			 (in_seq->data[i+2] == 'U' && in_seq->data[i+3] == 'A' && in_seq->data[i+4] == 'A') ||	 
			(in_seq->data[i+2] == 'T' && in_seq->data[i+3] == 'A' && in_seq->data[i+4] == 'G') ||
			 (in_seq->data[i+2] == 'U' && in_seq->data[i+3] == 'A' && in_seq->data[i+4] == 'G') ||
			(in_seq->data[i+2] == 'T' && in_seq->data[i+3] == 'G' && in_seq->data[i+4] == 'A') ||
			 (in_seq->data[i+2] == 'U' && in_seq->data[i+3] == 'G' && in_seq->data[i+4] == 'A') ) //Stop-Codon then end orf
				{
					orfs[4]->finish = 1; // end orf
					c=0; // wait for new ATG-Start-Codon
				}
				orfs[4]->sequence[(orfs[4]->length)] = translate(in_seq->data[i+2], in_seq->data[i+3], in_seq->data[i+4],&(orfs[4]->dna_num[(orfs[4]->length)]));	// translate in aminoacid and safe in orf
				++orfs[4]->length; 
			}
			else 
			{
				if((in_seq->data[i+2] == 'T' && in_seq->data[i+3] == 'A' && in_seq->data[i+4] == 'A') ||
			 (in_seq->data[i+2] == 'U' && in_seq->data[i+3] == 'A' && in_seq->data[i+4] == 'A') ||
			(in_seq->data[i+2] == 'T' && in_seq->data[i+3] == 'A' && in_seq->data[i+4] == 'G') ||
			 (in_seq->data[i+2] == 'U' && in_seq->data[i+3] == 'A' && in_seq->data[i+4] == 'G') ||
			(in_seq->data[i+2] == 'T' && in_seq->data[i+3] == 'G' && in_seq->data[i+4] == 'A') ||
			 (in_seq->data[i+2] == 'U' && in_seq->data[i+3] == 'G' && in_seq->data[i+4] == 'A') )
				{
					orfs[5]->finish = 1;
					c=0;
				}
				orfs[5]->sequence[(orfs[5]->length)] = translate(in_seq->data[i+2], in_seq->data[i+3], in_seq->data[i+4],&(orfs[5]->dna_num[(orfs[5]->length)]) );	
				++orfs[5]->length;
			}		
		}

/****************************************************************************************/
/****************************************************************************************/	

		if (d == 1)
		{
			if(!(orfs[6]->finish))
				
									/* Da Crick-Strand :    TAA -> ATT
															UAA -> AUU
															TAG -> ATC
															UAG -> AUC
															TGA -> ACT
															UGA -> ACU */

			{
				if((in_seq->data[help-(i+2)] == 'A' && in_seq->data[help-(i+3)] == 'T' && in_seq->data[help-(i+4)] == 'T') ||
			 (in_seq->data[help-(i+2)] == 'A' && in_seq->data[help-(i+3)] == 'U' && in_seq->data[help-(i+4)] == 'U') ||	 
			(in_seq->data[help-(i+2)] == 'A' && in_seq->data[help-(i+3)] == 'T' && in_seq->data[help-(i+4)] == 'C') ||
			 (in_seq->data[help-(i+2)] == 'A' && in_seq->data[help-(i+3)] == 'U' && in_seq->data[help-(i+4)] == 'C') ||
			(in_seq->data[help-(i+2)] == 'A' && in_seq->data[help-(i+3)] == 'C' && in_seq->data[help-(i+4)] == 'T') ||
			 (in_seq->data[help-(i+2)] == 'A' && in_seq->data[help-(i+3)] == 'C' && in_seq->data[help-(i+4)] == 'U') ) //Stop-Codon then end orf
				{
					orfs[6]->finish = 1; // end orf
					d=0; // wait for new ATG-Start-Codon
				}
				orfs[6]->sequence[(orfs[6]->length)] = translate(inverse(in_seq->data[help-(i+2)]), inverse(in_seq->data[help-(i+3)]), inverse(in_seq->data[help-(i+4)]),&(orfs[6]->dna_num[(orfs[6]->length)]) );	// translate in aminoacid and safe in orf
				++orfs[6]->length; 
			}
			else 
			{
				if((in_seq->data[help-(i+2)] == 'A' && in_seq->data[help-(i+3)] == 'T' && in_seq->data[help-(i+4)] == 'T') ||
			 (in_seq->data[help-(i+2)] == 'A' && in_seq->data[help-(i+3)] == 'U' && in_seq->data[help-(i+4)] == 'U') ||
			(in_seq->data[help-(i+2)] == 'A' && in_seq->data[help-(i+3)] == 'T' && in_seq->data[help-(i+4)] == 'C') ||
			 (in_seq->data[help-(i+2)] == 'A' && in_seq->data[help-(i+3)] == 'U' && in_seq->data[help-(i+4)] == 'C') ||
			(in_seq->data[help-(i+2)] == 'A' && in_seq->data[help-(i+3)] == 'C' && in_seq->data[help-(i+4)] == 'T') ||
			 (in_seq->data[help-(i+2)] == 'A' && in_seq->data[help-(i+3)] == 'C' && in_seq->data[help-(i+4)] == 'U') )
				{
					orfs[7]->finish = 1;
					d=0;
				}
				orfs[7]->sequence[(orfs[7]->length)] = translate(inverse(in_seq->data[help-(i+2)]), inverse(in_seq->data[help-(i+3)]), inverse(in_seq->data[help-(i+4)]),&(orfs[7]->dna_num[(orfs[7]->length)]) );
				++orfs[7]->length;
			}		
		}

/****************************************************************************************/	

		if (e == 1)
		{
			if(!(orfs[8]->finish))

			{
				if((in_seq->data[help-(i+1)] == 'A' && in_seq->data[help-(i+2)] == 'T' && in_seq->data[help-(i+3)] == 'T') ||
			 (in_seq->data[help-(i+1)] == 'A' && in_seq->data[help-(i+2)] == 'U' && in_seq->data[help-(i+3)] == 'U') ||	 
			(in_seq->data[help-(i+1)] == 'A' && in_seq->data[help-(i+2)] == 'T' && in_seq->data[help-(i+3)] == 'C') ||
			 (in_seq->data[help-(i+1)] == 'A' && in_seq->data[help-(i+2)] == 'U' && in_seq->data[help-(i+3)] == 'C') ||
			(in_seq->data[help-(i+1)] == 'A' && in_seq->data[help-(i+2)] == 'C' && in_seq->data[help-(i+3)] == 'T') ||
			 (in_seq->data[help-(i+1)] == 'A' && in_seq->data[help-(i+2)] == 'C' && in_seq->data[help-(i+3)] == 'U') ) //Stop-Codon then end orf
				{
					orfs[8]->finish = 1; // end orf
					e=0; // wait for new ATG-Start-Codon
				}
				orfs[8]->sequence[(orfs[8]->length)] = translate(inverse(in_seq->data[help-(i+1)]), inverse(in_seq->data[help-(i+2)]), inverse(in_seq->data[help-(i+3)]),&(orfs[8]->dna_num[(orfs[8]->length)]) );	// translate in aminoacid and safe in orf
				++orfs[8]->length; 
			}
			else 
			{
				if((in_seq->data[help-(i+1)] == 'A' && in_seq->data[help-(i+2)] == 'T' && in_seq->data[help-(i+3)] == 'T') ||
			 (in_seq->data[help-(i+1)] == 'A' && in_seq->data[help-(i+2)] == 'U' && in_seq->data[help-(i+3)] == 'U') ||
			(in_seq->data[help-(i+1)] == 'A' && in_seq->data[help-(i+2)] == 'T' && in_seq->data[help-(i+3)] == 'C') ||
			 (in_seq->data[help-(i+1)] == 'A' && in_seq->data[help-(i+2)] == 'U' && in_seq->data[help-(i+3)] == 'C') ||
			(in_seq->data[help-(i+1)] == 'A' && in_seq->data[help-(i+2)] == 'C' && in_seq->data[help-(i+3)] == 'T') ||
			 (in_seq->data[help-(i+1)] == 'A' && in_seq->data[help-(i+2)] == 'C' && in_seq->data[help-(i+3)] == 'U') )
				{
					orfs[9]->finish = 1;
					e=0;
				}
				orfs[9]->sequence[(orfs[9]->length)] = translate(inverse(in_seq->data[help-(i+1)]), inverse(in_seq->data[help-(i+2)]), inverse(in_seq->data[help-(i+3)]), &(orfs[9]->dna_num[(orfs[9]->length)]) );	
				++orfs[9]->length;
			}		
		}

/****************************************************************************************/	

		if (f == 1)
		{
			if(!(orfs[10]->finish))

			{
				if((in_seq->data[help-(i)] == 'A' && in_seq->data[help-(i+1)] == 'T' && in_seq->data[help-(i+2)] == 'T') ||
			 (in_seq->data[help-(i)] == 'A' && in_seq->data[help-(i+1)] == 'U' && in_seq->data[help-(i+2)] == 'U') ||	 
			(in_seq->data[help-(i)] == 'A' && in_seq->data[help-(i+1)] == 'T' && in_seq->data[help-(i+2)] == 'C') ||
			 (in_seq->data[help-(i)] == 'A' && in_seq->data[help-(i+1)] == 'U' && in_seq->data[help-(i+2)] == 'C') ||
			(in_seq->data[help-(i)] == 'A' && in_seq->data[help-(i+1)] == 'C' && in_seq->data[help-(i+2)] == 'T') ||
			 (in_seq->data[help-(i)] == 'A' && in_seq->data[help-(i+1)] == 'C' && in_seq->data[help-(i+2)] == 'U') ) //Stop-Codon then end orf
				{
					orfs[10]->finish = 1; // end orf
					f=0; // wait for new ATG-Start-Codon
				}
				orfs[10]->sequence[(orfs[10]->length)] = translate(inverse(in_seq->data[help-(i)]), inverse(in_seq->data[help-(i+1)]), inverse(in_seq->data[help-(i+2)]), &(orfs[10]->dna_num[(orfs[10]->length)]));	// translate in aminoacid and safe in orf
				++orfs[10]->length; 
			}
			else 
			{
				if((in_seq->data[help-(i)] == 'A' && in_seq->data[help-(i+1)] == 'T' && in_seq->data[help-(i+2)] == 'T') ||
			 (in_seq->data[help-(i)] == 'A' && in_seq->data[help-(i+1)] == 'U' && in_seq->data[help-(i+2)] == 'T') ||
			(in_seq->data[help-(i)] == 'A' && in_seq->data[help-(i+1)] == 'T' && in_seq->data[help-(i+2)] == 'C') ||
			 (in_seq->data[help-(i)] == 'A' && in_seq->data[help-(i+1)] == 'U' && in_seq->data[help-(i+2)] == 'C') ||
			(in_seq->data[help-(i)] == 'A' && in_seq->data[help-(i+1)] == 'C' && in_seq->data[help-(i+2)] == 'T') ||
			 (in_seq->data[help-(i)] == 'A' && in_seq->data[help-(i+1)] == 'C' && in_seq->data[help-(i+2)] == 'U') )
				{
					orfs[11]->finish = 1;
					f=0;
				}
				orfs[11]->sequence[(orfs[11]->length)] = translate(inverse(in_seq->data[help-(i)]), inverse(in_seq->data[help-(i+1)]), inverse(in_seq->data[help-(i+2)]), &(orfs[11]->dna_num[(orfs[11]->length)]));	
				++orfs[11]->length;
			}		
		}

/****************************************************************************************/
/****************************************************************************************/
		if(((in_seq->data[i] == 'A' && in_seq->data[i+1] == 'T' && in_seq->data[i+2] == 'G') && (a==0)) ||
			( (in_seq->data[i] == 'A' && in_seq->data[i+1] == 'U' && in_seq->data[i+2] == 'G') && (a==0)))
		{
			a=1; // start found in first reading frame
			flag_no_start=0;
			
			if(!(orfs[0]->finish)) // orf not yet started in orfs[0]
			{
				orfs[0]->length=0;  
				orfs[0]->sequence[(orfs[0]->length)] = translate(in_seq->data[i],in_seq->data[i+1],in_seq->data[i+2], &(orfs[0]->dna_num[(orfs[0]->length)]) );
				++orfs[0]->length;
			}

			else if(!(orfs[1]->finish))// orf not yet started in orfs[1]
			{
				orfs[1]->length=0;
				orfs[1]->sequence[(orfs[1]->length)] = translate(in_seq->data[i],in_seq->data[i+1],in_seq->data[i+2], &(orfs[1]->dna_num[(orfs[1]->length)]) );
				++orfs[1]->length;
			}
			else  // orf finished in orfs[1] and orfs[0], so find the max and delete the min
			{
				if((orfs[0]->length) >= (orfs[1]->length)) // orfs[0] = Max
				{
					free(orfs[1]->sequence);
					free(orfs[1]->dna_num);
					free(orfs[1]);
					if((orfs[1] = calloc(1, sizeof(struct orf))) == NULL) // free memory
						error("orf_finder(): Out of memory!");
					if((orfs[1]->sequence = calloc(j+1, sizeof(char))) == NULL)
						error("orf_finder(): Out of memory!");
					if((orfs[1]->dna_num = calloc(j+1, sizeof(char))) == NULL)
						error("orf_finder(): Out of memory!");
					orfs[1]->length=0;							// new init
					orfs[1]->finish = 0;
 					orfs[1]->sequence[(orfs[1]->length)] = translate(in_seq->data[i],in_seq->data[i+1],in_seq->data[i+2], &(orfs[1]->dna_num[(orfs[1]->length)]) );
					++orfs[1]->length;
				}
				else
				{
					free(orfs[0]->dna_num);
					free(orfs[0]->sequence);
					free(orfs[0]);
					if((orfs[0] = calloc(1, sizeof(struct orf))) == NULL) // free memory
						error("orf_finder(): Out of memory!");
					if((orfs[0]->sequence = calloc(j+1, sizeof(char))) == NULL)
						error("orf_finder(): Out of memory!");
					if((orfs[0]->dna_num = calloc(j+1, sizeof(char))) == NULL)
						error("orf_finder(): Out of memory!");
					orfs[0]->length = 0;						// new init
					orfs[0]->finish = 0;
					orfs[0]->sequence[(orfs[0]->length)] = translate(in_seq->data[i],in_seq->data[i+1],in_seq->data[i+2], &(orfs[0]->dna_num[(orfs[0]->length)]) );
					++orfs[0]->length;
				}	
			}
		}

/****************************************************************************************/
		
		if(((in_seq->data[i+1] == 'A' && in_seq->data[i+2] == 'T' && in_seq->data[i+3] == 'G') && (b==0)) ||
			( (in_seq->data[i+1] == 'A' && in_seq->data[i+2] == 'U' && in_seq->data[i+3] == 'G') && (b==0)))
		{
			b=1; // start found in second reading frame
			flag_no_start=0;

			if(!(orfs[2]->finish)) // orf not yet started in orfs[2]
			{
				orfs[2]->length=0;                          
				orfs[2]->sequence[(orfs[2]->length)] = translate(in_seq->data[i+1],in_seq->data[i+2],in_seq->data[i+3], &(orfs[2]->dna_num[(orfs[2]->length)]) );
				++orfs[2]->length;
			}

			else if(!(orfs[3]->finish))// orf not yet started in orfs[3]
			{
				orfs[3]->length=0;
				orfs[3]->sequence[(orfs[3]->length)] = translate(in_seq->data[i+1],in_seq->data[i+2],in_seq->data[i+3], &(orfs[3]->dna_num[(orfs[3]->length)]) );
				++orfs[3]->length;
			}
			else  // orf finished in orfs[3] and orfs[2], so find the max and delete the min
			{
				if((orfs[2]->length) >= (orfs[3]->length)) // orfs[2] = Max
				{
					free(orfs[3]->dna_num);
					free(orfs[3]->sequence);
					free(orfs[3]);
					if((orfs[3] = calloc(1, sizeof(struct orf))) == NULL) // free memory
						error("orf_finder(): Out of memory!");
					if((orfs[3]->sequence = calloc(j+1, sizeof(char))) == NULL)
						error("orf_finder(): Out of memory!");
					if((orfs[3]->dna_num = calloc(j+1, sizeof(char))) == NULL)
						error("orf_finder(): Out of memory!");

					orfs[3]->length=0;							// new init
					orfs[3]->finish = 0;
 					orfs[3]->sequence[(orfs[3]->length)] = translate(in_seq->data[i+1],in_seq->data[i+2],in_seq->data[i+3], &(orfs[3]->dna_num[(orfs[3]->length)]) );
					++orfs[3]->length;
				}
				else
				{
					free(orfs[2]->dna_num);
					free(orfs[2]->sequence);
					free(orfs[2]);
					if((orfs[2] = calloc(1, sizeof(struct orf))) == NULL) // free memory
						error("orf_finder(): Out of memory!");
					if((orfs[2]->sequence = calloc(j+1, sizeof(char))) == NULL)
						error("orf_finder(): Out of memory!");	
					if((orfs[2]->dna_num = calloc(j+1, sizeof(char))) == NULL)
						error("orf_finder(): Out of memory!");

					orfs[2]->length = 0;						// new init
					orfs[2]->finish = 0;
					orfs[2]->sequence[(orfs[2]->length)] = translate(in_seq->data[i+1],in_seq->data[i+2],in_seq->data[i+3], &(orfs[2]->dna_num[(orfs[2]->length)]) );
					++orfs[2]->length;
				}	
			}
		}
/****************************************************************************************/
		
		if(((in_seq->data[i+2] == 'A' && in_seq->data[i+3] == 'T' && in_seq->data[i+4] == 'G') && (c==0)) ||
			( (in_seq->data[i+2] == 'A' && in_seq->data[i+3] == 'U' && in_seq->data[i+4] == 'G') && (c==0) ))
		{
			c=1; // start found in third reading frame
			flag_no_start=0;

			if(!(orfs[4]->finish)) // orf not yet started in orfs[4]
			{
				orfs[4]->length=0;                          
				orfs[4]->sequence[(orfs[4]->length)] = translate(in_seq->data[i+2],in_seq->data[i+3],in_seq->data[i+4], &(orfs[4]->dna_num[(orfs[4]->length)]) );
				++orfs[4]->length;
			}

			else if(!(orfs[5]->finish))// orf not yet started in orfs[5]
			{
				orfs[5]->length=0;
				orfs[5]->sequence[(orfs[5]->length)] = translate(in_seq->data[i+2],in_seq->data[i+3],in_seq->data[i+4], &(orfs[5]->dna_num[(orfs[5]->length)]) );
				++orfs[5]->length;
			}
			else  // orf finished in orfs[5] and orfs[4], so find the max and delete the min
			{
				if((orfs[4]->length) >= (orfs[5]->length)) // orfs[4] = Max
				{
					free(orfs[5]->dna_num);
					free(orfs[5]->sequence);
					free(orfs[5]);
					if((orfs[5] = calloc(1, sizeof(struct orf))) == NULL) // free memory
						error("orf_finder(): Out of memory!");
					if((orfs[5]->sequence = calloc(j+1, sizeof(char))) == NULL)
						error("orf_finder(): Out of memory!");
					if((orfs[5]->dna_num = calloc(j+1, sizeof(char))) == NULL)
						error("orf_finder(): Out of memory!");
					orfs[5]->length=0;							// new init
					orfs[5]->finish = 0;
 					orfs[5]->sequence[(orfs[5]->length)] = translate(in_seq->data[i+2],in_seq->data[i+3],in_seq->data[i+4], &(orfs[5]->dna_num[(orfs[5]->length)]) );
					++orfs[5]->length;
				}
				else
				{
					free(orfs[4]->dna_num);
					free(orfs[4]->sequence);
					free(orfs[4]);
					if((orfs[4] = calloc(1, sizeof(struct orf))) == NULL) // free memory
						error("orf_finder(): Out of memory!");
					if((orfs[4]->sequence = calloc(j+1, sizeof(char))) == NULL)
						error("orf_finder(): Out of memory!");
					if((orfs[4]->dna_num = calloc(j+1, sizeof(char))) == NULL)
						error("orf_finder(): Out of memory!");
					orfs[4]->length = 0;						// new init
					orfs[4]->finish = 0;
					orfs[4]->sequence[(orfs[4]->length)] = translate(in_seq->data[i+2],in_seq->data[i+3],in_seq->data[i+4], &(orfs[4]->dna_num[(orfs[4]->length)]) );
					++orfs[4]->length;
				}	
			}
		}

/****************************************************************************************/
/****************************************************************************************/


		if(((in_seq->data[help-(i+2)] == 'T' && in_seq->data[help-(i+3)] == 'A' && in_seq->data[help-(i+4)] == 'C') && (d==0)) ||
			( (in_seq->data[help-(i+2)] == 'U' && in_seq->data[help-(i+3)] == 'A' && in_seq->data[help-(i+4)] == 'C') && (d==0) ))
		{
			d=1; // start found in last but two reading frame
			flag_no_start=0;

			if(!(orfs[6]->finish)) // orf not yet started in orfs[6]
			{
				orfs[6]->length=0;                          
				orfs[6]->sequence[(orfs[6]->length)] = translate(in_seq->data[help-(i+2)],in_seq->data[help-(i+3)],in_seq->data[help-(i+4)], &(orfs[6]->dna_num[(orfs[6]->length)]) );
				++orfs[6]->length;
			}

			else if(!(orfs[7]->finish))// orf not yet started in orfs[7]
			{
				orfs[7]->length=0;
				orfs[7]->sequence[(orfs[7]->length)] = translate(in_seq->data[help-(i+2)],in_seq->data[help-(i+3)],in_seq->data[help-(i+4)], &(orfs[7]->dna_num[(orfs[7]->length)]) );
				++orfs[7]->length;
			}
			else  // orf finished in orfs[7] and orfs[6], so find the max and delete the min
			{
				if((orfs[6]->length) >= (orfs[7]->length)) // orfs[6] = Max
				{
					free(orfs[7]->dna_num);
					free(orfs[7]->sequence);
					free(orfs[7]);
					if((orfs[7] = calloc(1, sizeof(struct orf))) == NULL) // free memory
						error("orf_finder(): Out of memory!");
					if((orfs[7]->sequence = calloc(j+1, sizeof(char))) == NULL)
						error("orf_finder(): Out of memory!");
					if((orfs[7]->dna_num = calloc(j+1, sizeof(char))) == NULL)
						error("orf_finder(): Out of memory!");
					orfs[7]->length=0;							// new init
					orfs[7]->finish = 0;
 					orfs[7]->sequence[(orfs[7]->length)] = translate(in_seq->data[help-(i+2)],in_seq->data[help-(i+3)],in_seq->data[help-(i+4)], &(orfs[7]->dna_num[(orfs[7]->length)]) );
					++orfs[7]->length;
				}
				else
				{
					free(orfs[6]->dna_num);
					free(orfs[6]->sequence);
					free(orfs[6]);
					if((orfs[6] = calloc(1, sizeof(struct orf))) == NULL) // free memory
						error("orf_finder(): Out of memory!");
					if((orfs[6]->sequence = calloc(j+1, sizeof(char))) == NULL)
						error("orf_finder(): Out of memory!");
					if((orfs[6]->dna_num = calloc(j+1, sizeof(char))) == NULL)
						error("orf_finder(): Out of memory!");
					orfs[6]->length = 0;						// new init
					orfs[6]->finish = 0;
					orfs[6]->sequence[(orfs[6]->length)] = translate(in_seq->data[help-(i+2)],in_seq->data[help-(i+3)],in_seq->data[help-(i+4)], &(orfs[6]->dna_num[(orfs[6]->length)]) );
					++orfs[6]->length;
				}	
			}
		}

/****************************************************************************************/

		if(((in_seq->data[help-(i+1)] == 'T' && in_seq->data[help-(i+2)] == 'A' && in_seq->data[help-(i+3)] == 'C') && (e==0)) ||
			( (in_seq->data[help-(i+1)] == 'U' && in_seq->data[help-(i+2)] == 'A' && in_seq->data[help-(i+3)] == 'C') && (e==0)))
		{
			e=1; // start found in last but one reading frame
			flag_no_start=0;

			if(!(orfs[8]->finish)) // orf not yet started in orfs[8]
			{
				orfs[8]->length=0;                          
				orfs[8]->sequence[(orfs[8]->length)] = translate(in_seq->data[help-(i+1)],in_seq->data[help-(i+2)],in_seq->data[help-(i+3)], &(orfs[8]->dna_num[(orfs[8]->length)]) );
				++orfs[8]->length;
			}

			else if(!(orfs[9]->finish))// orf not yet started in orfs[9]
			{
				orfs[9]->length=0;
				orfs[9]->sequence[(orfs[9]->length)] = translate(in_seq->data[help-(i+1)],in_seq->data[help-(i+2)],in_seq->data[help-(i+3)], &(orfs[9]->dna_num[(orfs[9]->length)]) );
				++orfs[9]->length;
			}
			else  // orf finished in orfs[9] and orfs[8], so find the max and delete the min
			{
				if((orfs[8]->length) >= (orfs[9]->length)) // orfs[8] = Max
				{
					free(orfs[9]->dna_num);
					free(orfs[9]->sequence);
					free(orfs[9]);
					if((orfs[9] = calloc(1, sizeof(struct orf))) == NULL) // free memory
						error("orf_finder(): Out of memory!");
					if((orfs[9]->sequence = calloc(j+1, sizeof(char))) == NULL)
						error("orf_finder(): Out of memory!");
					if((orfs[9]->dna_num = calloc(j+1, sizeof(char))) == NULL)
						error("orf_finder(): Out of memory!");
					orfs[9]->length=0;							// new init
					orfs[9]->finish = 0;
 					orfs[9]->sequence[(orfs[9]->length)] = translate(in_seq->data[help-(i+1)],in_seq->data[help-(i+2)],in_seq->data[help-(i+3)], &(orfs[9]->dna_num[(orfs[9]->length)]) );
					++orfs[9]->length;
				}
				else
				{
//					printf(" im else also orf[0]<orf[1]\n"); // orfs[8] = Min
					free(orfs[8]->dna_num);
					free(orfs[8]->sequence);
					free(orfs[8]);
					if((orfs[8] = calloc(1, sizeof(struct orf))) == NULL) // free memory
						error("orf_finder(): Out of memory!");
					if((orfs[8]->sequence = calloc(j+1, sizeof(char))) == NULL)
						error("orf_finder(): Out of memory!");
					if((orfs[8]->dna_num = calloc(j+1, sizeof(char))) == NULL)
						error("orf_finder(): Out of memory!");
					orfs[8]->length = 0;						// new init
					orfs[8]->finish = 0;
					orfs[8]->sequence[(orfs[8]->length)] = translate(in_seq->data[help-(i+1)],in_seq->data[help-(i+2)],in_seq->data[help-(i+3)], &(orfs[8]->dna_num[(orfs[8]->length)]) );
					++orfs[8]->length;
				}	
			}
		}

/****************************************************************************************/

		if(((in_seq->data[help-i] == 'T' && in_seq->data[help-(i+1)] == 'A' && in_seq->data[help-(i+2)] == 'C') && (f==0)) ||
			( (in_seq->data[help-i] == 'U' && in_seq->data[help-(i+1)] == 'A' && in_seq->data[help-(i+2)] == 'C') && (f==0)))
		{
			f=1; // start found in last reading frame
			flag_no_start=0;

			if(!(orfs[10]->finish)) // orf not yet started in orfs[10]
			{
				orfs[10]->length=0;                          
				orfs[10]->sequence[(orfs[10]->length)] = translate(in_seq->data[help-(i)],in_seq->data[help-(i+1)],in_seq->data[help-(i+2)], &(orfs[10]->dna_num[(orfs[10]->length)]) );
				++orfs[10]->length;
			}

			else if(!(orfs[11]->finish))// orf not yet started in orfs[11]
			{
				orfs[11]->length=0;
				orfs[11]->sequence[(orfs[11]->length)] = translate(in_seq->data[help-(i)],in_seq->data[help-(i+1)],in_seq->data[help-(i+2)], &(orfs[11]->dna_num[(orfs[11]->length)]) );
				++orfs[11]->length;
			}
			else  // orf finished in orfs[11] and orfs[10], so find the max and delete the min
			{
				if((orfs[10]->length) >= (orfs[11]->length)) // orfs[10] = Max
				{
					free(orfs[11]->dna_num);
					free(orfs[11]->sequence);
					free(orfs[11]);
					if((orfs[11] = calloc(1, sizeof(struct orf))) == NULL) // free memory
						error("orf_finder(): Out of memory!");
					if((orfs[11]->sequence = calloc(j+1, sizeof(char))) == NULL)
						error("orf_finder(): Out of memory!");
					if((orfs[11]->dna_num = calloc(j+1, sizeof(char))) == NULL)
						error("orf_finder(): Out of memory!");
					orfs[11]->length=0;							// new init
					orfs[11]->finish = 0;
 					orfs[11]->sequence[(orfs[11]->length)] = translate(in_seq->data[help-(i)],in_seq->data[help-(i+1)],in_seq->data[help-(i+2)], &(orfs[11]->dna_num[(orfs[11]->length)]) );
					++orfs[11]->length;
				}
				else
				{
					free(orfs[10]->dna_num);
					free(orfs[10]->sequence);
					free(orfs[10]);
					if((orfs[10] = calloc(1, sizeof(struct orf))) == NULL) // free memory
						error("orf_finder(): Out of memory!");
					if((orfs[10]->sequence = calloc(j+1, sizeof(char))) == NULL)
						error("orf_finder(): Out of memory!");
					if((orfs[10]->dna_num = calloc(j+1, sizeof(char))) == NULL)
						error("orf_finder(): Out of memory!");
					orfs[10]->length = 0;						// new init
					orfs[10]->finish = 0;
					orfs[10]->sequence[(orfs[10]->length)] = translate(in_seq->data[help-(i)],in_seq->data[help-(i+1)],in_seq->data[help-(i+2)], &(orfs[10]->dna_num[(orfs[10]->length)]) );
					++orfs[10]->length;
				}	
			}
		}

/****************************************************************************************/
/****************************************************************************************/

		i = i+2;
	}
	// if last codon in reading frame is no stop codon -> translate last codon if available
	if (((help+1) % 3) == 0) // length = last index + 1 ;translate last codon for orf 1 and orf 6
	{
//		a:
		if (!orfs[0]->finish && orfs[0]->length) // open orf && already started before
		{
			orfs[0]->sequence[(orfs[0]->length)] = translate(in_seq->data[help-2], in_seq->data[help-1], in_seq->data[help], &(orfs[0]->dna_num[(orfs[0]->length)]));
			++orfs[0]->length;
		}
		else if (!orfs[1]->finish && orfs[1]->length)
		{
			orfs[1]->sequence[(orfs[1]->length)] = translate(in_seq->data[help-2], in_seq->data[help-1], in_seq->data[help], &(orfs[1]->dna_num[(orfs[1]->length)]));
			++orfs[1]->length;
		}
		
//		f:
		if (!(orfs[10]->finish) && orfs[10]->length )
		{
			orfs[10]->sequence[(orfs[10]->length)] = translate(inverse(in_seq->data[2]), inverse(in_seq->data[1]), inverse(in_seq->data[0]), &(orfs[10]->dna_num[(orfs[10]->length)]) );
			++orfs[10]->length;
		}
		else if (!orfs[11]->finish && orfs[11]->length )
		{
			orfs[11]->sequence[(orfs[11]->length)] = translate(inverse(in_seq->data[2]), inverse(in_seq->data[1]), inverse(in_seq->data[0]), &(orfs[11]->dna_num[(orfs[11]->length)]) );
			++orfs[11]->length;
		}
	}


	else if (((help+1) % 3) == 1)
	{
//		a:
		if (!orfs[0]->finish && orfs[0]->length) // open orf && already started before
		{
			orfs[0]->sequence[(orfs[0]->length)] = translate(in_seq->data[help-3], in_seq->data[help-2], in_seq->data[help-1], &(orfs[0]->dna_num[(orfs[0]->length)]) );
			++orfs[0]->length;
		}
		else if (!orfs[1]->finish && orfs[1]->length)
		{
			orfs[1]->sequence[(orfs[1]->length)] = translate(in_seq->data[help-3], in_seq->data[help-2], in_seq->data[help-1], &(orfs[1]->dna_num[(orfs[1]->length)]) );
			++orfs[1]->length;
		}
		

//		b:
		if (!orfs[2]->finish && orfs[2]->length) // open orf && already started before
		{
			orfs[2]->sequence[(orfs[2]->length)] = translate(in_seq->data[help-2], in_seq->data[help-1], in_seq->data[help], &(orfs[2]->dna_num[(orfs[2]->length)]) );
			++orfs[2]->length;
		}
		else if (!orfs[3]->finish && orfs[3]->length)
		{
			orfs[3]->sequence[(orfs[3]->length)] = translate(in_seq->data[help-2], in_seq->data[help-1], in_seq->data[help], &(orfs[3]->dna_num[(orfs[3]->length)]) );
			++orfs[3]->length;
		}
	

//		e:
		if (!(orfs[8]->finish) && orfs[8]->length )
		{
			orfs[8]->sequence[(orfs[8]->length)] = translate(inverse(in_seq->data[2]), inverse(in_seq->data[1]), inverse(in_seq->data[0]), &(orfs[8]->dna_num[(orfs[8]->length)]) );
			++orfs[8]->length;
		}
		else if (!orfs[9]->finish && orfs[9]->length )
		{
			orfs[9]->sequence[(orfs[9]->length)] = translate(inverse(in_seq->data[2]), inverse(in_seq->data[1]), inverse(in_seq->data[0]), &(orfs[9]->dna_num[(orfs[9]->length)]) );
			++orfs[9]->length;
		}
		

//		f:
		if (!(orfs[10]->finish) && orfs[10]->length )
		{
			orfs[10]->sequence[(orfs[10]->length)] = translate(inverse(in_seq->data[3]), inverse(in_seq->data[2]), inverse(in_seq->data[1]), &(orfs[10]->dna_num[(orfs[10]->length)]) );
			++orfs[10]->length;
		}
		else if (!orfs[11]->finish && orfs[11]->length )
		{
			orfs[11]->sequence[(orfs[11]->length)] = translate(inverse(in_seq->data[3]), inverse(in_seq->data[2]), inverse(in_seq->data[1]), &(orfs[11]->dna_num[(orfs[11]->length)]) );
			++orfs[11]->length;
		}
	}
	else
	{
//		a:
		if (!orfs[0]->finish && orfs[0]->length) // open orf && already started before
		{
			orfs[0]->sequence[(orfs[0]->length)] = translate(in_seq->data[help-4], in_seq->data[help-3], in_seq->data[help-2], &(orfs[0]->dna_num[(orfs[0]->length)]) );
			++orfs[0]->length;
		}
		else if (!orfs[1]->finish && orfs[1]->length)
		{
			orfs[1]->sequence[(orfs[1]->length)] = translate(in_seq->data[help-4], in_seq->data[help-3], in_seq->data[help-2], &(orfs[1]->dna_num[(orfs[1]->length)]));
			++orfs[1]->length;
		}
		

//		b:
		if (!orfs[2]->finish && orfs[2]->length) // open orf && already started before
		{
			orfs[2]->sequence[(orfs[2]->length)] = translate(in_seq->data[help-3], in_seq->data[help-2], in_seq->data[help-1], &(orfs[2]->dna_num[(orfs[2]->length)]) );
			++orfs[2]->length;
		}
		else if (!orfs[3]->finish && orfs[3]->length)
		{
			orfs[3]->sequence[(orfs[3]->length)] = translate(in_seq->data[help-3], in_seq->data[help-2], in_seq->data[help-1], &(orfs[3]->dna_num[(orfs[3]->length)]) );
			++orfs[3]->length;
		}
		

//		c:
		if (!orfs[4]->finish && orfs[4]->length) // open orf && already started before
		{
			orfs[4]->sequence[(orfs[4]->length)] = translate(in_seq->data[help-2], in_seq->data[help-1], in_seq->data[help],&(orfs[4]->dna_num[(orfs[4]->length)]) );
			++orfs[4]->length;
		}
		else if (!orfs[5]->finish && orfs[5]->length)
		{
			orfs[5]->sequence[(orfs[5]->length)] = translate(in_seq->data[help-2], in_seq->data[help-1], in_seq->data[help],&(orfs[5]->dna_num[(orfs[5]->length)]) );
			++orfs[5]->length;
		}
		

//		d:
		if (!(orfs[6]->finish) && orfs[6]->length )
		{
			orfs[6]->sequence[(orfs[6]->length)] = translate(inverse(in_seq->data[2]), inverse(in_seq->data[1]), inverse(in_seq->data[0]), &(orfs[6]->dna_num[(orfs[6]->length)]) );
			++orfs[6]->length;
		}
		else if (!orfs[7]->finish && orfs[7]->length )
		{
			orfs[7]->sequence[(orfs[7]->length)] = translate(inverse(in_seq->data[2]), inverse(in_seq->data[1]), inverse(in_seq->data[0]), &(orfs[7]->dna_num[(orfs[7]->length)]) );
			++orfs[7]->length;
		}
		

//		e:
		if (!(orfs[8]->finish) && orfs[8]->length )
		{
			orfs[8]->sequence[(orfs[8]->length)] = translate(inverse(in_seq->data[3]), inverse(in_seq->data[2]), inverse(in_seq->data[1]), &(orfs[8]->dna_num[(orfs[8]->length)]) );
			++orfs[8]->length;
		}
		else if (!orfs[9]->finish && orfs[9]->length )
		{
			orfs[9]->sequence[(orfs[9]->length)] = translate(inverse(in_seq->data[3]), inverse(in_seq->data[2]), inverse(in_seq->data[1]), &(orfs[9]->dna_num[(orfs[9]->length)]) );
			++orfs[9]->length;
		}
		

//		f:
		if (!(orfs[10]->finish) && orfs[10]->length )
		{
			orfs[10]->sequence[(orfs[10]->length)] = translate(inverse(in_seq->data[4]), inverse(in_seq->data[3]), inverse(in_seq->data[2]), &(orfs[10]->dna_num[(orfs[10]->length)]) );
			++orfs[10]->length;
		}
		else if (!orfs[11]->finish && orfs[11]->length )
		{
			orfs[11]->sequence[(orfs[11]->length)] = translate(inverse(in_seq->data[4]), inverse(in_seq->data[3]), inverse(in_seq->data[2]), &(orfs[11]->dna_num[(orfs[11]->length)]) );
			++orfs[11]->length;
		}
		
	}

// now we got 12 orfs and we search for the longest
	if((ret_seq = calloc(1,sizeof(struct seq)))==NULL)
	{
		error("orf_finder(): Out of memory");
	}

	max_length=orfs[0]->length;
	max_orf=0;

	for (i=1; i != 12; ++i)
	{
		if(max_length<orfs[i]->length)
		{
			max_length=orfs[i]->length;
			max_orf=i;
		}
	}
	if((ret_seq->data = calloc(max_length+1, sizeof(char)))==NULL)
	{
		error("orf_finder(): Out of memory");
	}
	if((ret_seq->dna_num = calloc(max_length+1, sizeof(char)))==NULL)
	{
		error("orf_finder(): Out of memory");
	}
	if(max_orf >5)ret_seq->crick_strand=1;
	ret_seq->name=in_seq->name;
	ret_seq->length = max_length;
	ret_seq->num=in_seq->num;
	ret_seq->orf_frame=max_orf/2;     // As orf[1] && orf[0] = Reading frame 0 ,...
	for(i=0; i!= max_length; i++)
	{
		ret_seq->dna_num[i]=orfs[max_orf]->dna_num[i];
		ret_seq->data[i]=orfs[max_orf]->sequence[i];
	}

// free memory
	for (x = 0; x != 12; ++x)
	{
		free(orfs[x]->dna_num);
		free(orfs[x]->sequence);
		free(orfs[x]);
	}
	free(orfs);


			
	if(flag_no_start){ 
		char tmp[300];
		sprintf(tmp,"orf_finder(): No Startcodon found in Sequence %s !\n", in_seq->name );
		if(para->ORF_FRAME)
			printf("\n%sStarting translation for Sequence %s in reading frame 0\n\n", tmp,in_seq->name);
		else 
			error(tmp);
	}		
	return ret_seq;
}
