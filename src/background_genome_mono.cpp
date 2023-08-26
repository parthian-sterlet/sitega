#define _CRT_SECURE_NO_WARNINGS

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include  <time.h>
#define Min(a,b) ((a)>(b))? (b):(a);
#define Max(a,b) ((a)>(b))? (a):(b);
#define SEQLEN 12000
#define NCHR 50
#define NBIN 50

struct seqm {
	int num;
	int len;//length
	double fat;//length atgc
//	int nat;// % at	
	int don;// ready	
};
int compare_len(const void *X1, const void *X2)
{
	struct seqm *S1 = (struct seqm *)X1;
	struct seqm *S2 = (struct seqm *)X2;
	if (S1->len - S2->len > 0)return 1;
	if (S1->len - S2->len < 0)return -1;
	return 0;
}
int compare_fat(const void* X1, const void* X2)
{
	struct seqm* S1 = (struct seqm*)X1;
	struct seqm* S2 = (struct seqm*)X2;
	if (S1->fat - S2->fat > 0)return 1;
	if (S1->fat - S2->fat < 0)return -1;
	return 0;
}
int compare_num(const void* X1, const void* X2)
{
	struct seqm* S1 = (struct seqm*)X1;
	struct seqm* S2 = (struct seqm*)X2;
	if (S1->num - S2->num > 0)return 1;
	if (S1->num - S2->num < 0)return -1;
	return 0;
}

char *TransStr(char *d)
{
	int i, c, lens;
	lens = strlen(d);
	for (i = 0; i < lens; i++)
	{
		c = int(d[i]);
		if (c < 97) d[i] = char(c + 32);
		//else break;
	}
	return(d);
}
void DelChar(char *str, char c)
{
	int i, lens, size;

	size = 0;
	lens = strlen(str);
	for (i = 0; i < lens; i++)
	{
		if (str[i] != c)str[size++] = str[i];
	}
	str[size] = '\0';
}
int CheckStr(char *file, char *d, int n, int print, int bad)
{
	int i, len, ret,di;
	len = strlen(d);
	ret = 1;	
	for (i = 0; i < len; i++)
	{
		di = (int)(d[i]);
		if (strchr("atgcATGC", di )!= NULL)continue;		
		if (strchr("nN", di) != NULL) bad++;
		else
		{
			printf("File %s; sequence %d position %d (%c) bad. Sequence too bad!\n", file, n, i + 1, d[i]);
			exit(1);
		}
	}
	if (bad>0)
	{
		//if (print == 1)printf("File %s; sequence %d, %d positions are bad. Sequence processed!\n", file, n, bad);
		return 1;
	}
	return(ret);
}
void EvalSeq(char *file, int &nseq, int olen)
{
	char l[SEQLEN], d[SEQLEN], head[400];
	int fl = 0;
	FILE  *in;

	if ((in = fopen(file, "rt")) == NULL)
	{
		printf("Input file %s can't be opened!\n", file);
		exit(1);
	}
	char symbol = fgetc(in);
	rewind(in);
	nseq = 0;
	int n = 0;
	memset(d, 0, sizeof(d));
	while (n >= 0)
	{
		if (fgets(l, sizeof(l), in) == NULL) fl = -1;
		if (*l == '\n' && fl != -1)continue;
		if (((*l == symbol) || (fl == -1)) && (fl != 0))
		{
			int lenx = (int)strlen(d);
			int bad = 0;
			int check = CheckStr(file, d, n, 1,bad);
			lenx -= bad;
			if (lenx >= olen && check == 1)nseq++;
			n++;
			if (fl == -1)
			{
				fclose(in);
				break;
			}
		}
		if (*l == symbol)
		{
			memset(head, 0, sizeof(head));
			DelChar(l, '\n');
			DelChar(l, '\r');
			strcpy(head, l);
			fl = 0; continue;
		}
		if (fl == 0)
		{
			memset(d, 0, sizeof(d));
			DelChar(l, '\n');
			DelChar(l, '\r');
			strcpy(d, l);
			fl = 1; continue;
		}
		if (strlen(d) + strlen(l) > sizeof(d))
		{
			printf("Size is large...");
			printf("l:%s\nstrlen(l):%zu\n", l, strlen(l));
			printf("d:%s\nstrlen(d):%zu\n", d, strlen(d));
			exit(1);
		}
		DelChar(l, '\n');
		DelChar(l, '\r');
		strcat(d, l);
	}
}
void EvalLen(char *file, int *len, int *bad, int olen)
{
	char l[SEQLEN], d[SEQLEN], head[400];
	int fl = 0;
	FILE  *in;

	if ((in = fopen(file, "rt")) == NULL)
	{
		printf("Input file %s can't be opened!\n", file);
		exit(1);
	}
	char symbol = fgetc(in);
	rewind(in);
	memset(d, 0, sizeof(d));
	int nn = 0, n = 0;
	while (n >= 0)
	{
		if (fgets(l, sizeof(l), in) == NULL) fl = -1;
		if (*l == '\n' && fl != -1)continue;
		if (((*l == symbol) || (fl == -1)) && (fl != 0))
		{
			int lenx = (int)strlen(d);
			int bad1 = 0;
			int check = CheckStr(file, d, n, 0, bad1);
			lenx -= bad1;
			if (lenx >= olen && check == 1)
			{
				len[n] = lenx;
				bad[n] = bad1;
				n++;
			}
			nn++;
			if (fl == -1)
			{
				fclose(in);
				break;
			}
		}
		if (*l == symbol)
		{
			memset(head, 0, sizeof(head));
			DelChar(l, '\n');
			DelChar(l, '\r');
			strcpy(head, l);
			fl = 0; continue;
		}
		if (fl == 0)
		{
			memset(d, 0, sizeof(d));
			DelChar(l, '\n');
			DelChar(l, '\r');
			strcpy(d, l);
			fl = 1; continue;
		}
		if (strlen(d) + strlen(l) > sizeof(d))
		{
			printf("Size is large...");
			printf("l:%s\nstrlen(l):%zu\n", l, strlen(l));
			printf("d:%s\nstrlen(d):%zu\n", d, strlen(d));
			exit(1);
		}
		DelChar(l, '\n');
		DelChar(l, '\r');
		strcat(d, l);
	}
}
char *TransStrBack(char *d)
{
	int i, c, lens;
	lens = strlen(d);
	for (i = 0; i < lens; i++)
	{
		c = int(d[i]);
		if (c >= 97) d[i] = char(c - 32);
		//else break;
	}
	return(d);
}
int ComplStr(char *d)
{
	int i, len;
	len = strlen(d);
	char d1[SEQLEN];
	strcpy(d1, d);
	//	memset(d,0,sizeof(d));
	for (i = 0; i < len; i++)
	{
		switch (d1[len - i - 1])
		{
		case 'A': {d[i] = 'T'; break; }
		case 'T': {d[i] = 'A'; break; }
		case 'C': {d[i] = 'G'; break; }
		case 'G': {d[i] = 'C'; break; }
		default: d[i] = 'N';
		}
	}
	return 1;
}
void ReadSeq(char *file, int nseq, int *len, char ***peak_real, int olen)
{
	char l[SEQLEN], d[2][SEQLEN], head[400];
	int fl = 0, j;
	FILE  *in;

	if ((in = fopen(file, "rt")) == NULL)
	{
		printf("Input file %s can't be opened!\n", file);
		exit(1);
	}
	char symbol = fgetc(in);
	rewind(in);
	int nn = 0, n = 0;
	while (n >= 0)
	{
		if (fgets(l, sizeof(l), in) == NULL) fl = -1;
		if (*l == '\n' && fl != -1)continue;
		if (((*l == symbol) || (fl == -1)) && (fl != 0))
		{
			int lenx = (int)strlen(d[0]);
			int bad = 0;
			int check = CheckStr(file, d[0], n, 0,bad);
			lenx -= bad;
			nn++;
			if (lenx >= olen && check == 1)
			{
				TransStrBack(d[0]);
				d[0][len[n]] = '\0';
				strcpy(d[1], d[0]);
				ComplStr(d[1]);
				d[1][len[n]] = '\0';
				for (j = 0; j < 2; j++)
				{
					strcpy(peak_real[j][n], d[j]);
					peak_real[j][n][len[n]] = '\0';
				}
				n++;
			}
			else
			{
			//	if (lenx < olen)printf("Short peak %d (Len %d) ignored\n", n + 1, lenx);
			//	if (check == -1)printf("Unusual symbol, peak %d ignored\n%s\n", n + 1, d[0]);
			}
			if (fl == -1)
			{
				fclose(in);
				break;
			}
		}
		if (*l == symbol)
		{
			memset(head, 0, sizeof(head));
			DelChar(l, '\n');
			DelChar(l, '\r');
			strcpy(head, l);
			fl = 0; continue;
		}
		if (fl == 0)
		{
			memset(d[0], 0, sizeof(d[0]));
			DelChar(l, '\n');
			DelChar(l, '\r');
			strcpy(d[0], l);
			fl = 1; continue;
		}
		if (strlen(d[0]) + strlen(l) > sizeof(d[0]))
		{
			printf("Size is large...");
			printf("l:%s\nstrlen(l):%zu\n", l, strlen(l));
			printf("d:%s\nstrlen(d):%zu\n", d[0], strlen(d[0]));
			exit(1);
		}
		DelChar(l, '\n');
		DelChar(l, '\r');
		strcat(d[0], l);
	}
}
int IdeLet(char c)
{
	int ret;
	switch (c) {
	case 'A': ret = 0; break;
	case 'C': ret = 1; break;
	case 'G': ret = 2; break;
	case 'T': ret = 3; break;
	case 'N': ret = -1; break;
	default: ret = -2;
	}
	return(ret);
}
void GetSost(char* d, int word, int size, int* sost)
{
	int i, j, k, i_sost, let;
	char letter[5] = "ACGT";
	int ten[8] = { 1, 4, 16, 64, 256, 1024, 4096, 16384 };
	int lens = (int)strlen(d);
	for (i = 0; i < size; i++)sost[i] = 0;
	int word1 = word - 1;
	int befo = 0, test = 0;
	for (i = 0; i < lens - word1; i++)
	{
		i_sost = 0;
		int coden = IdeLet(d[i + word1]);
		if (befo == 1 && coden >= 0)test = 1;
		if (befo == 0)
		{
			test = 1;
			for (j = 0; j < word; j++)
			{
				int codei = IdeLet(d[i + j]);
				if (codei < 0)
				{
					test = 0;
					break;
				}
			}
		}
		if (test == 1)
		{
			let = 0;
			for (j = word1; j >= 0; j--)
			{
				for (k = 0; k < 4; k++)
				{
					if (d[i + j] == letter[k]) { let = k; break; }
				}
				i_sost += ten[word - 1 - j] * let;
			}
			sost[i_sost]++;
		}
	}
}
int main(int argc, char *argv[])
{
	int i, j, k;
	char d[SEQLEN], d1[SEQLEN], filesta[10], fileend[10], genome[10];
	char filei[500], fileo1[500], fileosta_mo[500], fileosta_di[500], fileosta_mo_one[500], fileosta_di_one[500], file_log[500], filechr[NCHR][500], path_fasta[500];
	FILE *out, *in_seq[NCHR], *outm, *outd, *outm_one, *outd_one;
	if (argc != 14)
	{
		puts("Sintax: 1 path_genome 2file in_fa, 3file out_fa 4int height 5double mono prec 6int back_iter 7char genome (hg38 mm10 rn6 zf11 dm6 ce235 sc64 sch294 at10 gm21 zm73 mp61)");
		puts("8double stop_fraction 9file_out mono_for_vs_back 10 file_out di_for_back 11file_out mononucl statistics 12file_out dinucl statistics 13file_log");
		exit(1);
	}
	strcpy(filesta, "chr");
	strcpy(fileend, ".plain");
	char name_chr[NCHR][10];

	//human
	char name_chr_hg[24][3] = { "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y" };
	int n_chr_hg = 24;
	//hg19
	//	int sizelo[24]={249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,155270560,59373566};
	//hg38
	int sizelo_hg38[24] = { 248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468, 156040895, 57227415 };
	//mouse
	char name_chr_mm[21][3] = { "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "X", "Y" };
	int n_chr_mm = 21;	
	//mm9	
	//int sizelo[21]={197195432,181748087,159599783,155630120,152537259,149517037,152524553,131738871,124076172,129993255,121843856,121257530,120284312,125194864,103494974,98319150,95272651,90772031,61342430,166650296,15902555};	
	//mouse mm10
	int sizelo_mm10[21] = { 195471971, 182113224, 160039680, 156508116, 151834684, 149736546, 145441459, 129401213, 124595110, 130694993, 122082543, 120129022, 120421639, 124902244, 104043685, 98207768, 94987271, 90702639, 61431566, 171031299, 91744698 };
	//rat Rnor 6.0
	int n_chr_rn = 22;	
	int sizelo_rn6[22] = { 260522016, 249053267, 169034231, 182687754, 166875058, 140994061, 135012528, 123900184, 114175309, 107211142, 86241447, 46669029, 106807694, 104886043, 101769107, 84729064, 86533673, 83828827, 57337602, 54435887, 152453651, 18315841 };
	char name_chr_rn[22][3] = { "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "X", "Y" };
	//zebrafish 11
	char name_chr_zf[25][3] = {"1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25"};
	int n_chr_zf =25;
	//int sizelo_zf10[25]={58871917,59543403,62385949,76625712,71715914,60272633,74082188,54191831,56892771,45574255,45107271,49229541,51780250,51944548,47771147,55381981,53345113,51008593,48790377,55370968,45895719,39226288,46272358,42251103,36898761};
	int sizelo_zf11[25] = { 59578282, 59640629, 62628489, 78093715, 72500376, 60270059, 74282399, 54304671, 56459846, 45420867, 45484837, 49182954, 52186027, 52660232, 48040578, 55266484, 53461100, 51023478, 48449771, 55201332, 45934066, 39133080, 46223584, 42172926, 37502051 };
	//arabidopsis tair10
	char name_chr_at[5][3] = { "1", "2", "3", "4", "5" };
	int sizelo_at10[5] = { 30427671, 19698289, 23459830, 18585056, 26975502 };
	int n_chr_at = 5;
	//soybean Glycine_max_v2.1
	int n_chr_gm = 20;
	char name_chr_gm[20][3] = { "1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20" };
	int sizelo_gm21[20] = { 56831624, 48577505, 45779781, 52389146, 42234498, 51416486, 44630646, 47837940, 50189764, 51566898, 34766867, 40091314, 45874162, 49042192, 51756343, 37887014, 41641366, 58018742, 50746916, 47904181 };
	//Caenorhabditis elegans WBcel235
	char name_chr_ce[6][5] = { "I", "II", "III", "IV", "V", "X" };
	int sizelo_ce235[6] = { 15072434, 15279421, 13783801, 17493829, 20924180, 17718942 };
	int n_chr_ce = 6;
	//Drosophila melanogaster dm6
	char name_chr_dm[7][3] = { "2R", "2L", "3R", "3L", "X", "Y", "4" };	
	int sizelo_dm6[7] = { 25286936, 23513712, 32079331, 28110227, 23542271, 3667352, 1348131 };
	int n_chr_dm = 7;
	//yeast Schizosaccharomyces pombe ASM294v2 
	int sizelo_sch294[3] = { 5579133, 4539804, 2452883 };
	int n_chr_sch = 3;
	char name_chr_sch[3][5] = { "I", "II", "III" };
	//yeast Saccharomyces cerevisiae R64-1-1
	int sizelo_sc64[16] = { 230218, 813184, 316620, 1531933, 576874, 270161, 1090940, 562643, 439888, 745751, 666816, 1078177, 924431, 784333, 1091291, 948066 };
	int n_chr_sc = 16;
	char name_chr_sc[16][5] = { "I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI" };
	//Marchantia polymorpha 6.1
	int sizelo_mp61[10] = { 30584173, 29643427, 27142341, 26988051, 26794015, 23861560, 21963529, 21314552, 4523046, 7543715 };
	char name_chr_mp[10][3] = { "1","2","3","4","5","6","7","8","U","V" };
	int n_chr_mp = 10;
	//zea mays 73
	int sizelo_zm73[10]={301476924,237917468,232245527,242062272,217959525,169407836,176826311,175377492,157038028,149632204};
	char name_chr_zm[10][3] = {"1","2","3","4","5","6","7","8","9","10"};
	int n_chr_zm=10;
	int sizelo1[NCHR], sizelo2[NCHR];//lengths in bp, Mb
	int n_chr=0;
	strcpy(path_fasta, argv[1]);
	strcpy(filei, argv[2]);//in_file
	strcpy(fileo1, argv[3]);//out_file seq
	int height = atoi(argv[4]);
	double mono_prec = atof(argv[5]);
	int back_iter = atoi(argv[6]);
	strcpy(genome, argv[7]);
	double stop_thr = atof(argv[8]);//0.99;// fraction of peaks 100% covered with height background sequences
	strcpy(fileosta_mo, argv[9]);//out_file sta success rate
	strcpy(fileosta_di, argv[10]);//out_file sta success rate
	strcpy(fileosta_mo_one, argv[11]);//out_file sta dinucl content
	strcpy(fileosta_di_one, argv[12]);//out_file sta dinucl content
	strcpy(file_log, argv[13]);//out_file sta dinucl content	
	if (mono_prec >= 0.5 || mono_prec <= 0)
	{
		printf("Mononucleotide precision %f is wrong!\n", mono_prec);
		exit(1);
	}
	if (stop_thr > 1 || stop_thr <= 0)
	{
		printf("Fraction of peaks %f is wrong!\n", stop_thr);
		exit(1);
	}
	if (back_iter > 1000000 || back_iter <= 0)
	{
		printf("Number of iterations %d is wrong or too large!\n", back_iter);
		exit(1);
	}
	if (height > 1000 || height <= 0)
	{
		printf("Maximal number of background sequences per one foreground sequence: %d is wrong or too large!\n", height);
		exit(1);
	}
	int genome_rec = 0;
	if (strcmp(genome, "at10") == 0)
	{
		genome_rec = 1;
		n_chr = n_chr_at;
		for (i = 0; i < n_chr; i++)
		{
			sizelo1[i] = sizelo_at10[i];
			strcpy(name_chr[i], name_chr_at[i]);
		}
	}
	else
	{
		if (strcmp(genome, "hg38") == 0)
		{
			genome_rec = 1;
			n_chr = n_chr_hg;
			for (i = 0; i < n_chr; i++)
			{
				sizelo1[i] = sizelo_hg38[i];
				strcpy(name_chr[i], name_chr_hg[i]);
			}
		}
		else
		{
			if (strcmp(genome, "mm10") == 0)
			{
				genome_rec = 1;
				n_chr = n_chr_mm;
				for (i = 0; i < n_chr; i++)
				{
					sizelo1[i] = sizelo_mm10[i];
					strcpy(name_chr[i], name_chr_mm[i]);
				}
			}
			else
			{
				if (strcmp(genome, "rn6") == 0)
				{
					genome_rec = 1;
					n_chr = n_chr_rn;
					for (i = 0; i < n_chr; i++)
					{
						sizelo1[i] = sizelo_rn6[i];
						strcpy(name_chr[i], name_chr_rn[i]);
					}
				}
				else
				{
					if (strcmp(genome, "dm6") == 0)
					{
						genome_rec = 1;
						n_chr = n_chr_dm;
						for (i = 0; i < n_chr; i++)
						{
							sizelo1[i] = sizelo_dm6[i];
							strcpy(name_chr[i], name_chr_dm[i]);
						}
					}
					else
					{
						if (strcmp(genome, "ce235") == 0)
						{
							genome_rec = 1;
							n_chr = n_chr_ce;
							for (i = 0; i < n_chr; i++)
							{
								sizelo1[i] = sizelo_ce235[i];
								strcpy(name_chr[i], name_chr_ce[i]);
							}
						}
						else
						{
							if (strcmp(genome, "sc64") == 0)
							{
								genome_rec = 1;
								n_chr = n_chr_sc;
								for (i = 0; i < n_chr; i++)
								{
									sizelo1[i] = sizelo_sc64[i];
									strcpy(name_chr[i], name_chr_sc[i]);
								}
							}
							else
							{
								if (strcmp(genome, "sch294") == 0)
								{
									genome_rec = 1;
									n_chr = n_chr_sch;
									for (i = 0; i < n_chr; i++)
									{
										sizelo1[i] = sizelo_sch294[i];
										strcpy(name_chr[i], name_chr_sch[i]);
									}
								}
								else
								{
									if (strcmp(genome, "gm21") == 0)
									{
										genome_rec = 1;
										n_chr = n_chr_gm;
										for (i = 0; i < n_chr; i++)
										{
											sizelo1[i] = sizelo_gm21[i];
											strcpy(name_chr[i], name_chr_gm[i]);
										}
									}
									else
									{
										if (strcmp(genome, "mp61") == 0)
										{
											genome_rec = 1;
											n_chr = n_chr_mp;
											for (i = 0; i < n_chr; i++)
											{
												sizelo1[i] = sizelo_mp61[i];
												strcpy(name_chr[i], name_chr_mp[i]);
											}
										}
										else
										{
											if (strcmp(genome, "zm73") == 0)
											{
												genome_rec = 1;
												n_chr = n_chr_zm;
												for (i = 0; i < n_chr; i++)
												{
													sizelo1[i] = sizelo_zm73[i];
													strcpy(name_chr[i], name_chr_zm[i]);
												}
											}
											else
											{
												if (strcmp(genome, "zf11") == 0)
												{
													genome_rec = 1;
													n_chr = n_chr_zf;
													for (i = 0; i < n_chr; i++)
													{
														sizelo1[i] = sizelo_zf11[i];
														strcpy(name_chr[i], name_chr_zf[i]);
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	if (genome_rec == 0)
	{
		printf("Genome %s is not recognized\n",genome);
		exit(1);
	}
	int mil = 1000000;
	if (strcmp(genome, "sc64") == 0 || strcmp(genome, "sch294") == 0)mil = 100000;
	for (i = 0; i < n_chr; i++)sizelo2[i] = sizelo1[i] / mil;
	int win_gomol = 50;//to compare homology and min paek length	
	srand((unsigned)time(NULL));	
	{
		FILE* out_log;
		if ((out_log = fopen(file_log, "wt")) == NULL)
		{
			fprintf(out_log, "Error: Input file %s can't be opened!\n", file_log);
			exit(1);
		}
		fclose(out_log);
	}
	for (i = 0; i < n_chr; i++)
	{
		strcpy(filechr[i], path_fasta);
		strcat(filechr[i], filesta);
		strcat(filechr[i], name_chr[i]);
		strcat(filechr[i], fileend);
		if ((in_seq[i] = fopen(filechr[i], "rt")) == NULL)
		{
			printf("Input file %s can't be opened!\n", filechr[i]);
			exit(1);
		}
	}
	int tot_len = 0;
	for (i = 0; i < n_chr; i++)tot_len += sizelo2[i];
	int *len, nseq = 0, olen = win_gomol, *bad;
	EvalSeq(filei, nseq, olen);
	len = new int[nseq];
	if (len == NULL) { puts("Out of memory..."); exit(1); }
	bad = new int[nseq];
	if (bad == NULL){ puts("Out of memory..."); exit(1); }
	int dnseq = 0;	
	EvalLen(filei, len, bad, olen);
	int *hei;
	hei = new int[nseq];
	if (hei == NULL) { puts("Out of memory..."); exit(1); }
	char ***peak_real;
	peak_real = new char**[2];
	if (peak_real == NULL) { puts("Out of memory..."); exit(1); }
	for (i = 0; i < 2; i++)
	{
		peak_real[i] = new char*[nseq];
		for (j = 0; j < nseq; j++)
		{
			peak_real[i][j] = new char[len[j] + bad[j] + 1];
			if (peak_real[i][j] == NULL) { puts("Out of memory..."); exit(1); }
		}
	}
	ReadSeq(filei, nseq, len, peak_real, win_gomol);
	int heis = 0;
	{
		FILE* out_log;
		if ((out_log = fopen(file_log, "wt")) == NULL)
		{
			fprintf(out_log, "Input file %s can't be opened!\n", file_log);
			exit(1);
		}
		fprintf(out_log, "Required %d genomic sequences are found for %d input sequences out of total %d\n", height, heis, nseq);
		fclose(out_log);
	}
	int len_max, len_min;
	seqm *sort;
	sort = new seqm[nseq];
	if (sort == NULL) { puts("Out of memory..."); exit(1); }
	for (i = 0; i < nseq; i++)
	{		
		int mono_at = 0, mono_cg = 0;
		for (j = 0; j < len[i]; j++)
		{
			if ((peak_real[0][i][j] == 'A' || peak_real[0][i][j] == 'T') || (peak_real[0][i][j] == 'a' || peak_real[0][i][j] == 't'))mono_at++;
			else
			{
				if ((peak_real[0][i][j] == 'C' || peak_real[0][i][j] == 'G') || (peak_real[0][i][j] == 'c' || peak_real[0][i][j] == 'g'))mono_cg++;
			}
		}
		sort[i].num = i;
		sort[i].don = 0;
		sort[i].len = len[i];
	//	sort[i].nat = mono_at;				
		sort[i].fat = (double)mono_at/(mono_at + mono_cg);
	}
	for (i = 0; i < nseq; i++)hei[i] = 0;
	qsort((void*)(&sort[0]), nseq, sizeof(sort[0]), compare_len);
	len_max = sort[nseq - 1].len+1;
	len_min = sort[0].len;
	int pr_tot = 0;
	int fl = 0;
	int trys = nseq * back_iter / 5;
	int iter = 0;
	int height0 = 10;
	int nseqb = nseq * height0;
//	oliq *sele;
//	sele = new oliq[nseqb];
//	if (sele == NULL) { puts("Out of memory..."); exit(1); }
	int nseqb1 = nseqb - 1;
	int good = 0;
	int gomol = 0;
	/*
	while (iter < trys)
	{
		int rr = rand();
		int chr_z = 1, sum = 0, ra = rr % tot_len;
		for (i = 0; i < n_chr; i++)
		{
			sum += sizelo2[i];
			if (ra < sum)
			{
				chr_z = i;
				break;
			}
		}
		int z_len = sizelo2[chr_z];
		rr = rand();
		int rb = rr % z_len;
		rb *= mil;
		rr = rand();
		int rb1 = rr % 1000;
		rb1 *= 1000;
		rb += rb1;
		rr = rand();
		int rb2 = rr % 1000;
		rb += rb2;
		fseek(in_seq[chr_z], (long)(rb), SEEK_SET);
		int check = -1;
		iter++;
		while (fgets(d, len_max, in_seq[chr_z]) != NULL)
		{
			check = CheckStr(filechr[chr_z], d, 0, 0);
			if (check == 1)break;
		}
		TransStrBack(d);
		if (check == 1)
		{
			int gom = 0;
			for (i = 0; i < nseq; i++)
			{
				for (k = 0; k < 2; k++)
				{
					for (j = 0; j < len[i] - win_gomol + 1; j++)
					{
						if (strncmp(d, &peak_real[k][i][j], win_gomol) == 0)
						{
							gom = 1;
							break;
						}
					}
					if (gom == 1)break;
				}
				if (gom == 1)break;
			}
			if (gom == 1)
			{
				gomol++;
			}
			else
			{
				int cat = 0;
				for (i = 0; i < len_min; i++)
				{
					if ((d[i] == 'A' || d[i] == 'T') || (d[i] == 'a' || d[i] == 't')) cat++;
				}
				int len_cur = len_min;
				int done = 0;
				for (i = 0; i < nseq; i++)
				{
					if (sort[i].don >= height0)continue;
					for (j = len_cur; j < sort[i].len; j++)
					{
						if ((d[j] == 'A' || d[j] == 'T') || (d[j] == 'a' || d[j] == 't'))cat++;
					}
					double mono = fabs((double)(cat - sort[i].nat)) / sort[i].lena;
					if (mono < mono_prec)
					{
						strncpy(d1, d, sort[i].len);
						d1[sort[i].len] = '\0';
						sort[i].don++;
						if (sort[i].don == height0)good++;
						pr_tot++;
						done = 1;
						break;
					}
					len_cur = sort[i].len;
				}
			}
		}
		if (iter % 10000 == 0)
		{
			double rat = double(good) / nseq;
			printf("Iterations %5d\t Nseq_Background %5d\tFraction_Done %5f\tHomol %d\n", iter, pr_tot, rat, gomol);
			if (rat > stop_thr)break;
		}		
	}*/
	if ((out = fopen(fileo1, "wt")) == NULL)
	{
		printf("Input file %s can't be opened!", fileo1);
		exit(1);
	}	
	int size = 0;
	iter = pr_tot = 0;
	trys = nseq * back_iter;
	gomol = 0;
	int stop;
	{
		stop = (int)(stop_thr * nseq);
	}
	double step_fr = 1 / (double)NBIN;
	double fr_all_for[NBIN], fr_all_back[NBIN], fr_no[NBIN], val[NBIN];
	for (i = 0; i < NBIN; i++)fr_all_for[i] = fr_all_back[i] = fr_no[i] = 0;
	val[0] = step_fr;
	for (i = 1; i < NBIN; i++)val[i] = val[i - 1] + step_fr;
	int di[16], ditotback[16], ditotbak_len = 0;
	for (j = 0; j < 16; j++)ditotback[j] = di[j] = 0; 
	while (iter < trys && heis < nseq)
	{
		int rr = rand();
		int chr_z = 1, sum = 0, ra = rr % tot_len;
		for (i = 0; i < n_chr; i++)
		{
			sum += sizelo2[i];
			if (ra < sum)
			{
				chr_z = i;
				break;
			}
		}
		int z_len = sizelo2[chr_z];
		rr = rand();
		int rb = rr % z_len;
		rb *= 1000000;
		rr = rand();
		int rb1 = rr % 1000;
		rb1 *= 1000;
		rb += rb1;
		rr = rand();
		int rb2 = rr % 1000;
		rb += rb2;
		fseek(in_seq[chr_z], (long)(rb), SEEK_SET);
		int check = -1;
		iter++;
		while (fgets(d, len_max, in_seq[chr_z]) != NULL)
		{
			if (strstr(d, "N") != NULL || strstr(d, "n") != NULL)continue;
			else
			{
				check = 1;
				break;
			}
			//check = CheckStr(filechr[chr_z], d, 0, 0,len_max);
			//if (check == 1)break;
		}		
		if (check == 1)
		{
			TransStrBack(d);
			int gom = 0;
			for (i = 0; i < nseq; i++)
			{
				for (k = 0; k < 2; k++)
				{
					for (j = 0; j < len[i] - win_gomol + 1; j++)
					{
						if (strncmp(d, &peak_real[k][i][j], win_gomol) == 0)
						{
							gom = 1;
							break;
						}
					}
					if (gom == 1)break;
				}
				if (gom == 1)break;
			}
			if (gom == 1)
			{
				gomol++;
			}
			else
			{
				int cat = 0;
				for (i = 0; i < len_min; i++)if ((d[i] == 'A' || d[i] == 'T') || (d[i] == 'a' || d[i] == 't')) cat++;
				int len_cur = len_min;
				int done = 0;
				for (i = 0; i < nseq; i++)
				{
					if (hei[i] >= height)continue;
					for (j = len_cur; j < sort[i].len; j++)
					{
						if ((d[j] == 'A' || d[j] == 'T') || (d[j] == 'a' || d[j] == 't'))cat++;
					}
					double mono = fabs(sort[i].fat - (double)(cat) / sort[i].len);
					if (mono < mono_prec)
					{
						strncpy(d1, d, sort[i].len);
						d1[sort[i].len] = '\0';
						if (hei[i] < height)
						{
							hei[i]++;
							fprintf(out, ">peak%d_%d_Mo_%f\n", sort[i].num, hei[i], mono);
							fprintf(out, "%s\n", d1);
							if (hei[i] == height)heis++;
							done = 1;
							size++;
							double donj = (double)cat / sort[i].len;
							int jk = 0;
							for (k = 0; k < NBIN; k++)
							{
								if (donj <= val[k])
								{
									jk = k;
									break;
								}
							}
							fr_all_back[jk]++;
							ditotbak_len += sort[i].len - 1;
							GetSost(d1, 2, 16, di);
							for (j = 0; j < 16; j++)ditotback[j] += di[j];
						}						
					}
					len_cur = sort[i].len;
					if (done == 1)
					{
						pr_tot++;
						break;
					}
				}
			}
		}
		if (iter % 1000 == 0)
		{
			int inx = nseq - 1;
			for (i = nseq - 1; i >= 0; i--)
			{
				if (hei[i] < height)
				{
					len_max = sort[i].len+1;
					inx = i;
					break;
				}
			}
			//printf("Iterations %5d\t Nseq_Background %5d\tLenMax %d Inx %d Fraction_Done %5f\tHomol %d\n", iter, pr_tot, len_max, inx, (double)heis / nseq, gomol);
			if (iter % 10000 == 0)
			{
				FILE *out_log;
				if ((out_log = fopen(file_log, "wt")) == NULL)
				{
					fprintf(out_log, "Input file %s can't be opened!\n", file_log);
					exit(1);
				}
				fprintf(out_log, "Required %d genomic sequences are found for %d input sequences out of total %d\n", height, heis, nseq);
				fclose(out_log);
			}
			if (heis >= stop)break;			
		}
	}
	{
		FILE* out_log;
		if ((out_log = fopen(file_log, "wt")) == NULL)
		{
			fprintf(out_log, "Input file %s can't be opened!\n", file_log);
			exit(1);
		}
		fprintf(out_log, "Calculations are completed. Required %d genomes sequences are found for %d input sequences out of total %d\n", height, heis, nseq);
		fclose(out_log);
	}
	//qsort((void*)(&sele[0]), pr_tot, sizeof(sele[0]), compare_num);
	for (i = 0; i < n_chr; i++)fclose(in_seq[i]);
	fclose(out);
	for (i = 0; i < nseq; i++)sort[i].don = hei[i];
	delete[] hei;
	qsort((void*)(&sort[0]), nseq, sizeof(sort[0]), compare_num);		
	int success = 0, no_success = 0;
	if ((outm = fopen(fileosta_mo, "wt")) == NULL)
	{
		printf("Input file %s can't be opened!", fileosta_mo);
		exit(1);
	}
	if ((outd = fopen(fileosta_di, "wt")) == NULL)
	{
		printf("Input file %s can't be opened!", fileosta_di);
		exit(1);
	}
	if ((outm_one = fopen(fileosta_mo_one, "wt")) == NULL)
	{
		printf("Input file %s can't be opened!", fileosta_mo_one);
		exit(1);
	}
	if ((outd_one = fopen(fileosta_di_one, "wt")) == NULL)
	{
		printf("Input file %s can't be opened!", fileosta_di_one);
		exit(1);
	}
	for (i = 0; i < nseq; i++)
	{
		//fprintf(out1, "%d\t%d\t%f\n", sort[i].num + 1, sort[i].don, sort[i].fat);			
		int jk = 0;
		double donj = sort[i].fat;
		for (j = 0; j < NBIN; j++)
		{				
			if (donj <= val[j])
			{
				jk = j;
				break;
			}
		}
		if (sort[i].don >= height)success++;
		else
		{
			fr_no[jk]++;
			no_success++;
		}
		fr_all_for[jk]++;
	}
	for (i = 0; i < NBIN; i++)fr_all_for[i] /= nseq;		
	for (i = 0; i < NBIN; i++)fr_all_back[i] /= size;
	if (no_success > 0)
	{
		for (i = 0; i < NBIN; i++)fr_no[i] /= no_success;
	}
	/*for (i = 0; i < NBIN; i++)printf("\t%f", fr_no[i]);
	printf("\n");
	for (i = 0; i < NBIN; i++)printf("\t%f", fr_all_for[i]);
	printf("\n");*/
	int ista = 0, iend = NBIN - 1, dbin = NBIN / 10;
	for (i = 0; i < NBIN; i++)
	{
		int rest = (i + 1) % dbin;
		if (rest == 0)ista = i;
		if (fr_all_for[i] != 0 || fr_all_back[i] != 0)
		{
			break;
		}
	}
	for (i = NBIN - 1; i >= 0; i--)
	{
		int rest = (i + 1) % dbin;
		if (rest == 0)iend = i;
		if (fr_all_for[i] != 0 || fr_all_back[i] != 0)
		{
			break;
		}
	}
	char dinu[16][3] = { "AA","AC","AG","AT","CA","CC","CG","CT","GA","GC","GG","GT","TA","TC","TG","TT" };
	fprintf(outm_one, "A/T content (%%): selected foreground sequences vs. Foreground set\t%f\t%f\n", 100 * val[ista], 100 * val[iend]);
	fprintf(outm_one, "\tA/T content (%%)\n");
	fprintf(outd_one, "Dincucleotide frequencies (%%): selected foreground sequences vs. Foreground set\n");
	for (j = 0; j < 16; j++)fprintf(outd_one, "\t%s", dinu[j]);
	fprintf(outd_one, "\n");
	int ditot[16];
	int mo[4], motot[4];
	int count_tot = 0;
	double monotot = 0, lendtot=0, lenmtot=0;
	for (j = 0; j < 16; j++)di[j] = ditot[j]=0;
	for (j = 0; j < 4; j++)mo[j] = motot[j] = 0;
	for (i = 0; i < nseq; i++)
	{
		GetSost(peak_real[0][sort[i].num], 2, 16, di);
		GetSost(peak_real[0][sort[i].num], 1, 4, mo);
		int di_one = 0, mo_one = 0;
		for (j = 0; j < 16; j++)
		{
			di_one += di[j];
			ditot[j] += di[j];
		}		
		for (j = 0; j < 4; j++)
		{
			mo_one += mo[j];
			motot[j] += mo[j];
		}		
		monotot += mo[0]+mo[3];
		count_tot += sort[i].don;		
		lendtot += di_one;
		lenmtot += mo_one;
		if (sort[i].don < height)
		{					
			fprintf(outm_one, "#Seq %4d #FoundSeq %d", sort[i].num + 1, sort[i].don);					
			fprintf(outd_one, "#Seq %4d #FoundSeq %d", sort[i].num + 1, sort[i].don);
			for (j = 0; j < 16; j++)fprintf(outd_one, "\t%f", 100 * (double)di[j]/di_one);
			fprintf(outm_one, "\t%f", 100*(double)(mo[0] + mo[3]) / mo_one);
			fprintf(outd_one, "\n");
			fprintf(outm_one, "\n");
		}
	}
	monotot /= lenmtot;
	fprintf(outm_one, "AllSeq #AvFoundSeq %f", (double)count_tot/nseq);
	fprintf(outd_one, "AllSeq #AvFoundSeq %f", (double)count_tot / nseq);
	for (j = 0; j < 16; j++)fprintf(outd_one, "\t%f", 100*(double)ditot[j]/lendtot);
	
	fprintf(outm_one, "\t%f\n", 100 * monotot);
	fprintf(outd_one, "\n");									
	fclose(outm_one);
	fclose(outd_one);	
	{		
		fprintf(outm, "A/T content (%%): Foreground set vs. Background set\t%f\t%f\n\tForeground\tBackground\n", 100 * val[ista], 100 * val[iend]);
		for (i = ista; i <= iend; i++)
		{
			fprintf(outm, "%f\t%f\t%f\n", 100*val[i], 100*fr_all_for[i], 100*fr_all_back[i]);
		}
	}
	fclose(outm);
	fprintf(outd, "Dincucleotide frequencies (%%): Foreground set vs. Background set\n\tForeground\tBackground\n");	
	for (k = 0; k < 16; k++)
	{
		fprintf(outd, "%s\t%f\t%f\n", dinu[k], 100*(double)ditot[k] / lendtot, 100 * (double)ditotback[k] / ditotbak_len);
	}
	fclose(outd);
	for (k = 0; k < 2; k++)
	{
		for (i = 0; i < nseq; i++)
		{
			delete[] peak_real[k][i];
		}
		delete[] peak_real[k];
	}
	delete[] peak_real;
	delete[] len;
	delete[] bad;
	delete[] sort;
	return 0;
}