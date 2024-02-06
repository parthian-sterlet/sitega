#define _CRT_SECURE_NO_WARNINGS

#include  <stdio.h>
#include  <stdlib.h>
#include  <string.h>
#include  <math.h>
#include  <time.h>
#include  <ctype.h>

#define Min(a,b) ((a)>(b))? (b):(a);
#define Max(a,b) ((a)>(b))? (a):(b);
#define SEQLEN 5050
#define DIM 100
//double sost[DIM];
#define Min(a,b) ((a)>(b))? (b):(a);
#define Max(a,b) ((a)>(b))? (a):(b);
//double sost1[DIM];

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
char *TransStrBack(char *d)//a->A
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
int IdeLet(char c, char *alfabet)
{
	int i, ret = -1;
	for (i = 0; i < 4; i++)
	{
		if (c == alfabet[i]) { ret = i; break; }
	}
	return(ret);
}/*
void GetSost(char *d, int word, int *sost, char *letter)
{
int i, j, k, i_sost, let;
int ten[6]={1,4,16,64,256,1024};
int lens=strlen(d);
int size=1;
for(k=0;k<word;k++)size*=4;
for(i=0;i<size;i++)sost[i]=0;
for(i=0;i<lens-word+1;i++)
{
i_sost=0;
for(j=word-1;j>=0;j--)
{
for(k=0;k<4;k++)
{
if(d[i+j]==letter[k]){let=k;break;}
}
i_sost+=ten[word-1-j]*let;
}
sost[i_sost]++;
}
}
*/
int ComplStr(char *d)
{
	char *d1;
	int i, len;
	len = strlen(d);
	d1 = new char[len + 1];
	if (d1 == NULL) return 0;
	strcpy(d1, d);
	//	memset(d,0,sizeof(d));
	for (i = 0; i < len; i++)
	{
		switch (d1[len - i - 1])
		{
		case 'a': {d[i] = 't'; break; }
		case 't': {d[i] = 'a'; break; }
		case 'c': {d[i] = 'g'; break; }
		case 'g': {d[i] = 'c'; break; }
		case 'A': {d[i] = 'T'; break; }
		case 'T': {d[i] = 'A'; break; }
		case 'C': {d[i] = 'G'; break; }
		case 'G': {d[i] = 'C'; break; }
		case 'N': {d[i] = 'N'; break; }
		case 'n': {d[i] = 'n'; break; }
		default: d[i] = 'n';
		}
	}
	delete[] d1;
	return 1;
}
int CheckStr(char *d)
{
	int i, len, ret, size;
	ret = size = 0;
	len = strlen(d);
	for (i = 0; i < len; i++)
	{
		if (strchr("atgcn", (int)d[i]) != NULL) { continue; }
		else { ret++; }
	}
	return(ret);
}
int ConvertSym(int &c)
{
	char four[] = "atgc";
	char di[6][3] = { "ag", "tc", "at", "ac", "gt", "gc" };
	char tri[4][4] = { "agt", "agc", "gtc", "tca" };

	if ((c > 64) && (c < 91)) c = char(c + 32);

	switch (c)
	{
	case 'r': {
		c = (int)di[1][rand() % 2];
		break; }
	case 'y': {
		c = (int)di[1][rand() % 2];
		break; }
	case 'w': {
		c = (int)di[2][rand() % 2];
		break; }
	case 'm': {
		c = (int)di[3][rand() % 2];
		break; }
	case 'k': {
		c = (int)di[4][rand() % 2];
		break; }
	case 's': {
		c = (int)di[5][rand() % 2];
		break; }
	case 'd': {
		c = (int)tri[0][rand() % 3];
		break; }
	case 'v': {
		c = (int)tri[1][rand() % 3];
		break; }
	case 'b': {
		c = (int)tri[2][rand() % 3];
		break; }
	case 'h': {
		c = (int)tri[3][rand() % 3];
		break; }
	case 'n': {
		c = (int)four[rand() % 4];
		break; }
	default: {
		c = (int)'n';
		return -1; }
	}
	return 1;
}

int ReadSeq(char *file, int &n, int &len1, int &all_pos)
{
	char head[1000];
	int fl = 0, len;
	char symbol;
	int c;
	//	char cyfr[]="0123456789";
	FILE  *in;
	len1 = len = n = 0;
	if ((in = fopen(file, "rt")) == NULL)
	{
		printf("Input file %s can't be opened!\n", file);
		return -1;
	}
	symbol = fgetc(in);
	rewind(in);
	while ((c = fgetc(in)) != -1)
	{
		if ((char)c == symbol)
		{
			if (fgets(head, sizeof(head), in) == NULL)return -1;
			if (len > len1)len1 = len;
			len = 0;
			n++;
			continue;
		}
		if (strchr("\t\n ", c) != NULL)continue;
		if (strchr("ATGCatgc", c) != NULL)all_pos++;
		if (strchr("ATGCNatgcn", c) != NULL)
		{
			len++;
			continue;
		}
		else
		{
			printf("Unusual base: sequence N %d, letter position %d\n symbol %c\n%s", n, len, c, head);
			continue;
		}
	}
	if (len > len1)len1 = len;
	fclose(in);
	return 1;
}
void ReplaceChar(char *str, char c1, char c2)
{
	int i, len = strlen(str);
	for (i = 0; i < len; i++)
	{
		if (str[i] == c1) str[i] = c2;
	}
}
int StrNStr(char *str, char c, int n)
{
	int i, len = strlen(str);
	int k = 1;
	for (i = 0; i < len; i++)
	{
		if (str[i] == c)
		{
			if (k == n)return i;
			k++;
		}
	}
	return -1;
}
double UnderStol(char *str, int nstol, char razd)
{
	if (nstol == 0)return atof(str);
	char ret[100];
	memset(ret, 0, sizeof(ret));
	int p1 = StrNStr(str, razd, nstol);
	int p2 = StrNStr(str, razd, nstol + 1);
	if (p2 == -1)
	{
		p2 = strlen(str);
	}
	if (p1 == -1 || p2 == -1) return -1;
	int len = p2 - p1 - 1;
	strncpy(ret, &str[p1 + 1], len);
	ret[len] = '\0';
	return atof(ret);
}
int UnderStolStr(char *str, int nstol, char *ret, size_t size, char sep)
{
	memset(ret, '\0', size);
	int p1, p2, len;
	if (nstol == 0)
	{
		p2 = StrNStr(str, sep, 1);
		if (p2 == -1)p2 = strlen(str);
		strncpy(ret, str, p2);
		ret[p2] = '\0';
		return 1;
	}
	else
	{
		p1 = StrNStr(str, sep, nstol);
		p2 = StrNStr(str, sep, nstol + 1);
		if (p2 == -1)
		{
			p2 = strlen(str);
		}
		if (p1 == -1 || p2 == -1) return -1;
		len = p2 - p1 - 1;
		strncpy(ret, &str[p1 + 1], len);
		ret[len] = '\0';
		return 1;
	}
}
void Mix(double *a, double *b)
{
	double buf = *a;
	*a = *b;
	*b = buf;
}
struct due {
	double buf;
	int sta;
	int end;
	int num;
	void get_copy(due *a);
	//	void print_all(void);
};
void due::get_copy(due *a)
{
	a->num = num;
	a->sta = sta;
	a->buf = buf;
	a->end = end;
};
/*
void due::print_all(void)
{
printf("[%d;%d]%s\t", sta, end, s[num].oli);
}*/
//set of dinucleotides
struct city {
	char site[300];
	int size;
	int len;
	double c;
	double std;
	struct due tot[DIM];
	void get_copy(city *a);
	void sort_all(void);
	int get_file(char *file);
	//void city::fprint_tab(char *file);
}sta;
int city::get_file(char *file)
{
	FILE *in;
	if ((in = fopen(file, "rt")) == NULL)
	{
		printf("Input file %s can't be opened!", file);
		return -1;
	}
	char d[300];
	fgets(d, sizeof(d), in);
	DelChar(d,'\n');
	strcpy(site, d);
	fgets(d, sizeof(d), in);
	size = atoi(d);
	fgets(d, sizeof(d), in);
	len = atoi(d);
	fgets(d, sizeof(d), in);
	c = atof(d);
	std = 0.05;
	char sep = '\t', s[30];
	int i, test;
	for (i = 0; i < size; i++)
	{
		fgets(d, sizeof(d), in);
		tot[i].sta = atoi(d);
		test = UnderStolStr(d, 1, s, sizeof(s), sep);
		if (test == -1) { printf("Wrong format %s\n", d); return(-1); }
		tot[i].end = atoi(s);
		test = UnderStolStr(d, 2, s, sizeof(s), sep);
		if (test == -1) { printf("Wrong format %s\n", d); return(-1); }
		tot[i].buf = atof(s);
		test = UnderStolStr(d, 3, s, sizeof(s), sep);
		if (test == -1) { printf("Wrong format %s\n", d); return(-1); }
		tot[i].num = atoi(s);
	}
	fclose(in);
	return 1;
}
void city::get_copy(city *a)
{
	strcpy(a->site, site);
	a->size = size;
	a->std = std;
	a->len = len;
	a->c = c;
	int i;
	for (i = 0; i < size; i++)
	{
		tot[i].get_copy(&a->tot[i]);
	}
}
int compare_due(const void *X1, const void *X2)
{
	struct due *S1 = (struct due *)X1;
	struct due *S2 = (struct due *)X2;
	if (S1->sta - S2->sta > 0)return 1;
	if (S1->sta - S2->sta < 0)return -1;
	if (S1->end - S2->end > 0)return 1;
	if (S1->end - S2->end < 0)return -1;
	if (S1->num - S2->num > 0)return 1;
	if (S1->num - S2->num < 0)return -1;
	return 0;
}
void city::sort_all(void)
{
	qsort((void*)tot, size, sizeof(tot[0]), compare_due);
}

int main(int argc, char *argv[])
{
	int i, j, k, n;
	char head[1000], file_out_distt[300], file_out_distb[300], file_sitega[300], path_fasta[300], file_fasta[300];
	FILE *in, *out_distt, * out_distb;

	if (argc != 9)
	{
		printf("%s 1path_fasta 2sitega_matrix_file 3file_profile_fasta 4file out_dist_txt 5file out_dist_binary 6double pvalue_large 6double score_min 8double dpvalue", argv[0]);//5file out_cpp_arr 
		return -1;
	}
	char letter[] = "acgt";
	strcpy(path_fasta, argv[1]);
	strcpy(file_fasta, path_fasta);
	strcat(file_fasta, argv[3]);
	strcpy(file_sitega, argv[2]);
	strcpy(file_out_distt, argv[4]);
	strcpy(file_out_distb, argv[5]);
	double pvalue_large = atof(argv[6]);
	//strcpy(file_out_cpp_arr, argv[5]);
	double thr_bot = atof(argv[7]);
	double bin = atof(argv[8]);

	int nseq_pro = 0, len_pro = 0;
	int all_pos = 0;
	ReadSeq(file_fasta, nseq_pro, len_pro, all_pos);
	int nthr = 2 * (int)(pvalue_large*all_pos*1.05);
	double *thr;
	thr = new double[nthr];
	if (thr == NULL) { puts("Out of memory..."); return -1; }
	for (i = 0; i < nthr; i++)thr[i] = 0;
	int nthr_max = nthr - 1;
	char *dp;
	dp = new char[len_pro + 10];
	if (dp == NULL) { puts("Out of memory..."); return -1; }
	if ((in = fopen(file_fasta, "rt")) == NULL)
	{
		printf("Input file %s can't be opened!", file_fasta);
		return -1;
	}
	double all_pos_rec = 0;
	int count_val = 0;
	city sta;
	if(sta.get_file(file_sitega)==-1)
	{
		printf("Site %s function not found!", file_sitega);
		exit(1);
	}
	int len1 = sta.len;
	int rlen[DIM];
	for (j = 0; j < sta.size; j++)rlen[j] = (sta.tot[j].end - sta.tot[j].sta + 1);
	for (n = 0; n < nseq_pro; n++)
	{
		if (n % 100 == 0)printf("%5d %f\n", n, thr[nthr_max]);
		fgets(head, sizeof(head), in);
		memset(dp, 0, len_pro + 1);
		fgets(dp, len_pro + 2, in);
		DelChar(dp, '\n');
		TransStr(dp);
		int len_pro1 = strlen(dp);
		int gom = 1;
		for (i = 0; i < len_pro1; i++)
		{
			int di = (int)dp[i];
			if (strchr(letter, di) == NULL)
			{
				gom = 0;
				break;
			}
		}
		if (gom == 0)continue;
		int len21 = len_pro1 - len1;
		for (int compl1 = 0; compl1 < 2; compl1++)
		{
			if (compl1 == 1) if (ComplStr(dp) != 1) { puts("Out of memory..."); return -1; }
			char d2[SEQLEN];
			double p = -1000;
			for (i = 0; i <= len21; i++)
			{
				strncpy(d2, &dp[i], len1);
				d2[len1] = '\0';
				if (strstr(d2, "n") != NULL) { continue; }
				all_pos_rec++;
				double score = sta.c;
				for (j = 0; j < sta.size; j++)
				{
					double fm = 0;
					for (k = sta.tot[j].sta; k <= sta.tot[j].end; k++)
					{
						int cod = 4 * IdeLet(d2[k], letter) + IdeLet(d2[k + 1], letter);
						if (sta.tot[j].num == cod) { fm++; }
					}
					if (fm != 0)
					{
						fm /= rlen[j];
						score += sta.tot[j].buf*fm;
					}
				}
				score = 1 - fabs(score - 1);
				double thr_check = Max(thr_bot, thr[nthr_max]);
				if (score >= thr_check)
				{
					int gom = 0;
					for (j = 0; j < nthr; j++)
					{
						if (score >= thr[j])
						{
							//if (thr[j] != 0)
							{
								int ksta = Min(nthr_max, count_val);
								for (k = ksta; k > j; k--)
								{
									//									if (thr[k] == 0)continue;
									Mix(&thr[k - 1], &thr[k]);
								}
							}
							thr[j] = score;
							gom = 1;
							break;
						}
						if (gom == 1)break;
					}
					count_val++;
				}
			}
		}
	}
	fclose(in);
	int nthr_dist = 0;
	int nthr_final = nthr - 1;
	double fpr_pred = (double)1 / all_pos_rec;
	double thr_pred = thr[0];
	for (j = 1; j < nthr; j++)
	{
		double fpr = (double)(j + 1) / all_pos_rec;
		if ((thr[j] != thr_pred && fpr - fpr_pred > bin) || j == nthr_final)
		{
			nthr_dist++;
			if (fpr_pred >= pvalue_large)
			{
				break;
			}
			thr_pred = thr[j];
			fpr_pred = fpr;
		}
	}
	double *thr_dist, *fpr_dist;
	thr_dist = new double[nthr_dist];
	if (thr_dist == NULL) { puts("Out of memory..."); return -1; }
	fpr_dist = new double[nthr_dist];
	if (fpr_dist == NULL) { puts("Out of memory..."); return -1; }
	int count = 0;
	fpr_pred = (double)1 / all_pos_rec;
	thr_pred = thr[0];
	for (j = 1; j < nthr; j++)
	{
		double fpr = (double)(j + 1) / all_pos_rec;
		if ((thr[j] != thr_pred && fpr - fpr_pred > bin) || j == nthr_final)
		{
			if (j != nthr_final)
			{
				thr_dist[count] = thr_pred;
				fpr_dist[count] = fpr_pred;
			}
			else
			{
				thr_dist[count] = thr[j];
				fpr_dist[count] = fpr;
			}			
			fpr_dist[count] = -log10(fpr_dist[count]);
			count++;
			if (fpr_pred >= pvalue_large)
			{
				break;
			}
			thr_pred = thr[j];
			fpr_pred = fpr;
		}
	}	
	if ((out_distt = fopen(file_out_distt, "wt")) == NULL)
	{
		printf("Out file %s can't be opened!\n", file_out_distt);
		return -1;
	}
	for (j = 0; j < count; j++)fprintf(out_distt, "%.18f\t%.18g\n", thr_dist[j], fpr_dist[j]);
	fclose(out_distt);
	if ((out_distb = fopen(file_out_distb, "wb")) == NULL)
	{
		printf("Out file %s can't be opened!\n", file_out_distb);
		return -1;
	}
	fwrite(&sta, sizeof(sta), 1, out_distb);
	fwrite(&count, sizeof(int), 1, out_distb);
	fwrite(thr_dist, sizeof(double), count, out_distb);
	fwrite(fpr_dist, sizeof(double), count, out_distb);
	fclose(out_distb);
	FILE *out_sta;
	char file_sta[500];
	strcpy(file_sta,file_sitega); 
	strcat(file_sta,"_sta");
	if ((out_sta = fopen(file_sta, "wt")) == NULL)
	{
		printf("Out file %s can't be opened!\n", file_sta);
		return -1;
	}
	fprintf(out_sta, "%s\t%d\t", file_out_distt, nthr_dist);
	fprintf(out_sta, "%.18f\t%.18g\n", thr_dist[count-1], fpr_dist[count-1]);
	fclose(out_sta);
	delete[] thr_dist;
	delete[] fpr_dist;
	delete[] thr;
	delete[] dp;
	return 1;
}


