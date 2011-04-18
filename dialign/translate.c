/**
 *
 * translate.c:  
 *
 * 2004-08-30  Dorothea Emig Volker Menrad
 *             
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "struct.h"
#include "translate.h"
#include "parameters.h"
#include "orf.h"
#include "io.h"

void translate_sequence_collection_orf_frame(struct seq_col *in_seq_col)
{
	int i;
	struct seq_col* tmp;
	tmp = set_longest_orf(in_seq_col);
/**************
	
	reading frames :

	0: 123 123 123 123 123 ...
	1: X 123 123 123 123 123 ...
	2: XX 123 123 123 123 123 ...
	3: ... 321 321 321 321 321 XX
	4:  ... 321 321 321 321 321 X
	5:    ... 321 321 321 321 321 

**************/
	for(i=0; i!= in_seq_col->length; ++i){
		in_seq_col->seqs[i].orf_frame = tmp->seqs[i].orf_frame;
		translate_sequence( &(in_seq_col->seqs[i]) );
		free(tmp->seqs[i].data);
		free(tmp->seqs[i].dna_num);
	}
	
	free(tmp);
}


void translate_sequence_collection_default(struct seq_col *in_seq_col)
{
	int i;
	for(i=0; i!= in_seq_col->length; ++i)
		translate_sequence( &(in_seq_col->seqs[i]) );
}


void translate_sequence(struct seq *in_seq)
{
	int i, rest;
//	printf("%d\t%s\n%s\n%d\n\n\n",in_seq->num, in_seq->name, in_seq->data, in_seq->orf_frame );
	switch(in_seq->orf_frame)
	{
		case 0: break;
		case 1: --in_seq->length;
				in_seq->data=&(in_seq->data[1]);
				break;
		case 2: in_seq->length= in_seq->length-2;
				in_seq->data=&(in_seq->data[2]);
				break;
		case 3: in_seq->crick_strand=1;
				in_seq->length= in_seq->length-2;
				in_seq->data[in_seq->length]='\0';
				break;
		case 4: in_seq->crick_strand=1;
				--in_seq->length;
				in_seq->data[in_seq->length]='\0';
				break;
		case 5: in_seq->crick_strand=1;
				break;
	}

	struct seq* new;

	if(NULL == ( new = (calloc(1,sizeof(struct seq)))))
	{
		error("translate_sequence: Out of memory! ");
	}
	if(NULL == (new->data = (calloc(in_seq->length/3+1, sizeof(char)))))
	{
		error("translate_sequence: Out of memory! ");
	}
	if(NULL == (new->dna_num = (calloc(in_seq->length/3+1, sizeof(char)))))
	{
		error("translate_sequence: Out of memory! ");
	}

	rest = in_seq->length%3;

	for(i=0; i < (in_seq->length) - rest ; ++i)
	{
		int help = in_seq->length-1;
		
		if(in_seq->crick_strand==0)
		{
			new->data[new->length] = translate(in_seq->data[i], in_seq->data[i+1], in_seq->data[i+2], &(new->dna_num[new->length]) );
			++new->length;
		}
		else
		{
			new->data[new->length] = translate(inverse(in_seq->data[help - i]), inverse(in_seq->data[help - (i+1)]), inverse(in_seq->data[help - (i+2)]), &(new->dna_num[new->length]) );
			++new->length;
		}
		i=i+2;
	}
	
	switch(in_seq->orf_frame)
	{
		case 1:	--in_seq->data;
				free(in_seq->data);
				break;
		case 2: in_seq->data = in_seq->data-2;
				free(in_seq->data);
				break;
		default: free(in_seq->data);
	}

	in_seq->data = &(new->data[0]);
	in_seq->length = new->length;
	in_seq->dna_num = new->dna_num;

}

char translate(char first_base, char second_base, char third_base, char *dna_number)
{
	switch(first_base)
			{
				case 'A': 	switch(second_base)
							{
						 	case 'A':	switch(third_base)
										{
										case 'A':	*dna_number = 0;	// AAA
													return('K'); // Lysin = K	
										case 'T': 	*dna_number = 1;	// AAT
													return('N'); // Asparagin = N
										case 'U':	*dna_number = 2;	// AAU
													return('N'); // Asparagin = N
										case 'G': 	*dna_number  = 3;	// AAG
													return('K'); //Lysin = K
										case 'C': 	*dna_number = 4;	// AAC
													return('N'); //Asparagin = N
										default: 	error("No regular aminoacid !");
										}
							case 'T': 	switch(third_base)
										{
										case 'A':	*dna_number = 5;	// ATA
													return('I'); // Isoleucin = I
										case 'T': 	*dna_number = 6;	// ATT
													return('I'); // Isoleucin = I
										case 'G': 	*dna_number = 7;	// ATG
													return('M'); // Methionin = M
										case 'C': 	*dna_number = 8;	// ATC
													return('I'); // Isoleucin = I
										default: 	error("No regular aminoacid !");
										}
							case 'U':	switch(third_base)
										{
										case 'A':	*dna_number = 9;	// AUA
													return('I'); // Isoleucin = I	
										case 'U':	*dna_number = 10;	// AUU
													return('I'); // Isoleucin = I
										case 'G': 	*dna_number = 11;	// AUG
													return('M'); // Methionin= M
										case 'C': 	*dna_number = 12;	// AUC
													return('I'); // Isoleucin = I
										default: 	error("No regular aminoacid !");
										}
							case 'G': 	switch(third_base)
										{
										case 'A':	*dna_number = 13;	// AGA
													return('R'); // Arginin = R
										case 'T': 	*dna_number = 14;	// AGT
													return('S'); // Serin = S
										case 'U': 	*dna_number = 15;	// AGU
													return('S'); // Serin = S
										case 'G': 	*dna_number = 16;	// AGG
													return('R'); // Arginin = R
										case 'C': 	*dna_number = 17;	// AGC
													return('S'); // Serin = S
										default: 	error("No regular aminoacid !");
										}
							case 'C': 	switch(third_base)
										{
										case 'A': 	*dna_number = 18;	// ACA
													return('T'); // Threonin = T
										case 'T': 	*dna_number = 19;	// ACT
													return('T'); // Threonin = T
										case 'U':	*dna_number = 20;	// ACU
													return('T'); // Threonin = T
										case 'G':	*dna_number = 21;	// ACG
													return('T'); // Threonin = T
										case 'C':	*dna_number = 22;	// ACC
													return('T'); // Threonin = T
										default: 	error("No regular aminoacid !");
										}


							default: 	error("No regular aminoacid !");
							}
				case 'T': 	switch(second_base)
							{
						 	case 'A':	switch(third_base)
										{
										case 'A':	*dna_number = 23;	// TAA
													return('X'); // Stop
										case 'T': 	*dna_number = 24;	// TAT
													return('Y'); // Tyrosin
										case 'G': 	*dna_number = 25;	// TAG
													return('X'); // Stop
										case 'C': 	*dna_number = 26;	// TAC
													return('Y'); // Tyrosin
										default: 	error("No regular aminoacid !");
										}
							case 'T': 	switch(third_base)
										{
										case 'A':	*dna_number = 27;	// TTA
													return('L'); // Leucin
										case 'T': 	*dna_number = 28;	// TTT
													return('F'); // Phenylalanin
										case 'G': 	*dna_number = 29;	// TTG
													return('L'); // Leucin
										case 'C': 	*dna_number = 30;	// TTC
													return('F'); // Phenylalanin
										default: 	error("No regular aminoacid !");
										}
							case 'G': 	switch(third_base)
										{
										case 'A':	*dna_number = 31;	// TGA
													return('X'); // Stop
										case 'T': 	*dna_number = 32;	// TGT
													return('C'); // Cystein
										case 'G': 	*dna_number = 33;	// TGG
													return('W'); // Tryptophan
										case 'C': 	*dna_number = 34;	// TGC
													return('C'); // Cystein
										default: 	error("No regular aminoacid !");
										}
							case 'C': 	switch(third_base)
										{
										case 'A':	*dna_number = 35;	// TCA
													return('S'); // Serin
										case 'T':	*dna_number = 36;	// TCT
													return('S'); // Serin
										case 'G':	*dna_number = 37;	// TCG
													return('S'); // Serin
										case 'C':	*dna_number = 38;	// TCC
													return('S'); // Serin
										default: 	error("No regular aminoacid !");
										}
							default: 	error("no regular Aminoacid !");
							}
				case 'U': 	switch(second_base)
							{
						 	case 'A':	switch(third_base)
										{
										case 'A':	*dna_number = 39;	// UAA
													return('X'); // Stop
										case 'U': 	*dna_number = 40;	// UAU
													return('Y'); // Tyrosin
										case 'G': 	*dna_number = 41;	// UAG
													return('X'); // Stop
										case 'C': 	*dna_number = 42;	// UAC
													return('Y'); // Tyrosin
										default: 	error("No regular aminoacid !");
										}
							case 'U': 	switch(third_base)
										{
										case 'A':	*dna_number = 43;	// UUA
													return('L'); // Leucin
										case 'U': 	*dna_number = 44;	// UUU
													return('F'); // Phenylalanin
										case 'G': 	*dna_number = 45;	// UUG
													return('L'); // Leucin
										case 'C': 	*dna_number = 46;	// UUC
													return('F'); // Phenylalanin
										default: 	error("No regular aminoacid !");
										}
							case 'G': 	switch(third_base)
										{
										case 'A':	*dna_number = 47;	// UGA
													return('X'); // Stop
										case 'U': 	*dna_number = 48;	// UGU
													return('C'); // Cystein
										case 'G': 	*dna_number = 49;	// UGG
													return('W'); // Tryptophan
										case 'C': 	*dna_number = 50;	// UGC
													return('C'); // Cystein
										default: 	error("No regular aminoacid !");
										}
							case 'C': 	switch(third_base)
										{
										case 'A':	*dna_number = 51;	// UCA
													return('S'); // Serin
										case 'U':	*dna_number = 52;	// UCU
													return('S'); // Serin
										case 'G':	*dna_number = 53;	// UCG
													return('S'); // Serin
										case 'C':	*dna_number = 54;	// UCC
													return('S'); // Serin
										default: 	error("No regular aminoacid !");
										}
							default: 	error("no regular Aminoacid !");
							}
				case 'G': 	switch(second_base)
							{
						 	case 'A':	switch(third_base)
										{
										case 'A':	*dna_number = 55;	// GAA
													return('E'); // Glutaminsaeure
										case 'T': 	*dna_number = 56;	// GAT
													return('D'); // Asparaginsaeure
										case 'U': 	*dna_number = 57;	// GAU
													return('D'); // Asparaginsaeure
										case 'G': 	*dna_number = 58;	// GAG
													return('E'); // Glutaminsaeure
										case 'C': 	*dna_number = 59;	// GAC
													return('D'); // Asparaginsaeure
										default: 	error("No regular aminoacid !");
										}
							case 'T': 	switch(third_base)
										{
										case 'A':	*dna_number = 60;	// GTA
													return('V'); // Valin
										case 'T':	*dna_number = 61;	// GTT
													return('V'); // Valin
										case 'G':	*dna_number = 62;	// GTG
													return('V'); // Valin
										case 'C':	*dna_number = 63;	// GTC
													return('V'); // Valin
										default: 	error("No regular aminoacid !");
										}
							case 'U': 	switch(third_base)
										{
										case 'A':	*dna_number = 64;	// GUA
													return('V'); // Valin
										case 'U':	*dna_number = 65;	// GUU
													return('V'); // Valin
										case 'G':	*dna_number = 66;	// GUG
													return('V'); // Valin
										case 'C':	*dna_number = 67;	// GUC
													return('V'); // Valin
										default: 	error("No regular aminoacid !");
										}
							case 'G': 	switch(third_base)
										{
										case 'A':	*dna_number = 68;	// GGA
													return('G'); // Glycin
										case 'T':	*dna_number = 69;	// GGT
													return('G'); // Glycin
										case 'U':	*dna_number = 70;	// GGU
													return('G'); // Glycin
										case 'G':	*dna_number = 71;	// GGG
													return('G'); // Glycin
										case 'C':	*dna_number = 72;	// GGC
													return('G'); // Glycin
										default: 	error("No regular aminoacid !");
										}
							case 'C': 	switch(third_base)
										{
										case 'A':	*dna_number = 73;	// GCA
													return('A'); // Alanin
										case 'T':	*dna_number = 74;	// GCT
													return('A'); // Alanin
										case 'U':	*dna_number = 75;	// GCU
													return('A'); // Alanin
										case 'G':	*dna_number = 76;	// GCG
													return('A'); // Alanin
										case 'C':	*dna_number = 77;	// GCC
													return('A'); // Alanin
										default: 	error("No regular aminoacid !");
										}
							default: 	error("No regular aminoacid !");
							}
				case 'C': 	switch(second_base)
							{
						 	case 'A':	switch(third_base)
										{
										case 'A':	*dna_number = 78;	// CAA
													return('Q'); // Glutamin
										case 'T': 	*dna_number = 79;	// CAT
													return('H'); // Histidin
										case 'U': 	*dna_number = 80;	// CAU
													return('H'); // Histidin
										case 'G': 	*dna_number = 81;	// CAG
													return('Q'); // Glutamin
										case 'C': 	*dna_number = 82;	// CAC
													return('H'); // Histidin
										default: 	error("No regular aminoacid !");
										}
							case 'T': 	switch(third_base)
										{
										case 'A':	*dna_number = 83;	// CTA
													return('L'); // Leucin
										case 'T': 	*dna_number = 84;	// CTT
													return('L'); // Leucin
										case 'G': 	*dna_number = 85;	// CTG
													return('L'); // Leucin
										case 'C': 	*dna_number = 86;	// CTC
													return('L'); // Leucin
										default: 	error("No regular aminoacid !");
										}
							case 'U': 	switch(third_base)
										{
										case 'A':	*dna_number = 87;	// CUA
													return('L'); // Leucin
										case 'U': 	*dna_number = 88;	// CUU
													return('L'); // Leucin
										case 'G': 	*dna_number = 89;	// CUG
													return('L'); // Leucin
										case 'C': 	*dna_number = 90;	// CUC
													return('L'); // Leucin
										default: 	error("No regular aminoacid !");
										}
							case 'G': 	switch(third_base)
										{
										case 'A':	*dna_number = 91;	// CGA
													return('R'); // Arginin
										case 'T': 	*dna_number = 92;	// CGT
													return('R'); // Arginin
										case 'U': 	*dna_number = 93;	// CGU
													return('R'); // Arginin
										case 'G': 	*dna_number = 94;	// CGG
													return('R'); // Arginin
										case 'C': 	*dna_number = 95;	// CGC
													return('R'); // Arginin
										default: 	error("No regular aminoacid !");
										}
							case 'C': 		switch(third_base)
										{
										case 'A':	*dna_number = 96;	// CCA
													return('P'); // Prolin
										case 'T': 	*dna_number = 97;	// CCT
													return('P'); // Prolin
										case 'U': 	*dna_number = 98;	// CCU
													return('P'); // Prolin
										case 'G': 	*dna_number = 99;	// CCG
													return('P'); // Prolin
										case 'C': 	*dna_number = 100;	// CCC
													return('P'); // Prolin
										default: 	error("No regular aminoacid !");
										}
							default: 	error("No regular aminoacid !");
							}
				default: 	error("No regular aminoacid !");
			}
	abort();
}

char inverse(char base)
{
	switch(base)
	{
		case 'A':	return('T');
		case 'T':	return('A');
		case 'U':	return('A');
		case 'C':	return('G');
		case 'G':	return('C');
	}
	abort();
}

void retranslate_sequence(struct seq *in_seq)
{
	char *tmp_char= "000";
	int i;
	struct seq* tmp;
	if(NULL == ( tmp =(calloc(1,sizeof(struct seq)))))
	{
		error("retranslate_sequence(): Out of memory !");
	}
	if(NULL == ( tmp->data =(calloc((in_seq->length)*3+1,sizeof(char)))))
	{
		error("retranslate_sequence(): Out of memory !");
	}
	
	for(i=0; i!=in_seq->length; ++i)
	{
		tmp_char = retranslate(in_seq->dna_num[i]);
		strcat(tmp->data,tmp_char);
	}

	in_seq->length = strlen(tmp->data);
	free(in_seq->data);
	free(in_seq->dna_num);
	in_seq->data = tmp->data;
}

char* retranslate(char amino)
{
       	switch(amino)
	{
		case  0:	return "AAA";
		case  1:	return "AAT";
		case  2:	return "AAU";
		case  3:	return "AAG";
		case  4:	return "AAC";
		case  5:	return "ATA";
		case  6:	return "ATT";
		case  7:	return "ATG";
		case  8:	return "ATC";
		case  9:	return "AUA";
		case 10:	return "AUU";
		case 11:	return "AUG";
		case 12:	return "AUC";
		case 13:	return "AGA";
		case 14:	return "AGT";
		case 15:	return "AGU";
		case 16:	return "AGG";
		case 17:	return "AGC";
		case 18:	return "ACA";
		case 19:	return "ACT";
		case 20:	return "ACU";
		case 21:	return "ACG";
		case 22:	return "ACC";
		case 23:	return "TAA";
		case 24:	return "TAT";
		case 25:	return "TAG";
		case 26:	return "TAC";
		case 27:	return "TTA";
		case 28:	return "TTT";
		case 29:	return "TTG";
		case 30:	return "TTC";
		case 31:	return "TGA";
		case 32:	return "TGT";
		case 33:	return "TGG";
		case 34:	return "TGC";
		case 35:	return "TCA";
		case 36:	return "TCT";
		case 37:	return "TCG";
		case 38:	return "TCC";
		case 39:	return "UAA";
		case 40:	return "UAU";
		case 41:	return "UAG";
		case 42:	return "UAC";
		case 43:	return "UUA";
		case 44:	return "UUU";
		case 45:	return "UUG";
		case 46:	return "UUC";
		case 47:	return "UGA";
		case 48:	return "UGU";
		case 49:	return "UGG";
		case 50:	return "UGC";
		case 51:	return "UCA";
		case 52:	return "UCU";
		case 53:	return "UCG";
		case 54:	return "UCC";
		case 55:	return "GAA";
		case 56:	return "GAT";
		case 57:	return "GAU";
		case 58:	return "GAG";
		case 59:	return "GAC";
		case 60:	return "GTA";
		case 61:	return "GTT";
		case 62:	return "GTG";
		case 63:	return "GTC";
		case 64:	return "GUA";
		case 65:	return "GUU";
		case 66:	return "GUG";
		case 67:	return "GUC";
		case 68:	return "GGA";
		case 69:	return "GGT";
		case 70:	return "GGU";
		case 71:	return "GGG";
		case 72:	return "GGC";
		case 73:	return "GCA";
		case 74:	return "GCT";
		case 75:	return "GCU";
		case 76:	return "GCG";
		case 77:	return "GCC";
		case 78:	return "CAA";
		case 79:	return "CAT";
		case 80:	return "CAU";
		case 81:	return "CAG";
		case 82:	return "CAC";
		case 83:	return "CTA";
		case 84:	return "CTT";
		case 85:	return "CTG";
		case 86:	return "CTC";
		case 87:	return "CUA";
		case 88:	return "CUU";
		case 89:	return "CUG";
		case 90:	return "CUC";
		case 91:	return "CGA";
		case 92:	return "CGT";
		case 93:	return "CGU";
		case 94:	return "CGG";
		case 95:	return "CGC";
		case 96:	return "CCA";
		case 97:	return "CCT";
		case 98:	return "CCU";
		case 99:	return "CCG";
		case 100:	return "CCC";
		default: error("Sorry something wrong while retranslation !");
	}
	abort();
}
