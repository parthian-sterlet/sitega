#define _CRT_SECURE_NO_WARNINGS

#include  <stdio.h>
#include  <stdlib.h>
#include  <string.h>
#include  <ctype.h>
#include  <math.h>
#include  <time.h>
#define DIM 80
#define SEQLEN 30000
#define NCHR 50

int StrNStr(char* str, char c, int n)
{
	if (n == 0)return -1;
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
void DelHole(char* str)
{
	char* hole;
	hole = strstr(str, "\n");
	if (hole != NULL) *hole = 0;
}
int UnderStol(char* str, int nstol, char* ret, size_t size, char sep)
{
	memset(ret, 0, size);
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
struct ss {
	int num;
	char oli[3];
} s[16];
//#include "andy1.h"
#define TEN 10000
double sost[DIM];
#define Min(a,b) ((a)>(b))? (b):(a);
#define Max(a,b) ((a)>(b))? (a):(b);
double sost1[DIM];
double uw[DIM][DIM];
char** seq;
char** hseq;
int legr[DIM], regr[DIM], tssl[DIM], tssr[DIM];
char symbol;

struct due {
	double buf;
	int sta;
	int end;
	int num;
	void get_copy(due* a);
	void print_all(void);
};
void due::get_copy(due* a)
{
	a->num = num;
	a->sta = sta;
	a->buf = buf;
	a->end = end;
};
void due::print_all(void)
{
	printf("[%d;%d]%s\t", sta, end, s[num].oli);
}
//set of dinucleotides
struct city {
	char site[120];
	int size;
	int len;
	double min;
	double raz;
	struct due tot[DIM];
	void get_copy(city* a);
	int get_file(char* file);
	void sort_all(void);
	//void city::fprint_tab(char *file);
}sta;
int city::get_file(char* file)
{
	FILE* in;
	if ((in = fopen(file, "rt")) == NULL)
	{
		printf("Input file %s can't be opened!", file);
		return -1;
	}
	char d[120];
	fgets(d, sizeof(d), in);
	DelHole(d);
	strcpy(site, d);
	fgets(d, sizeof(d), in);
	size = atoi(d);
	fgets(d, sizeof(d), in);
	len = atoi(d);
	fgets(d, sizeof(d), in);
	min = atof(d);
	fgets(d, sizeof(d), in);
	raz = atof(d);
	char sep = '\t', s[30];
	int i, test;
	for (i = 0; i < size; i++)
	{
		fgets(d, sizeof(d), in);
		tot[i].sta = atoi(d);
		test = UnderStol(d, 1, s, sizeof(s), sep);
		if (test == -1) { printf("Wrong format %s\n", d); return(-1); }
		tot[i].end = atoi(s);
		test = UnderStol(d, 2, s, sizeof(s), sep);
		if (test == -1) { printf("Wrong format %s\n", d); return(-1); }
		tot[i].buf = atof(s);
		test = UnderStol(d, 3, s, sizeof(s), sep);
		if (test == -1) { printf("Wrong format %s\n", d); return(-1); }
		tot[i].num = atoi(s);
	}
	fclose(in);
	return 1;
}
void city::get_copy(city* a)
{
	strcpy(a->site, site);
	a->size = size;
	a->min = min;
	a->len = len;
	a->raz = raz;
	int i;
	for (i = 0; i < size; i++)
	{
		tot[i].get_copy(&a->tot[i]);
	}
}
int compare_due(const void* X1, const void* X2)
{
	struct due* S1 = (struct due*)X1;
	struct due* S2 = (struct due*)X2;
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
//#include "andy1_models.h"
//#include "andy1_get_model.h"

struct qbs {
	double q;
	double fn;
	int fp;
};
int compare_q(const void* X1, const void* X2)
{
	struct qbs* S1 = (struct qbs*)X1;
	struct qbs* S2 = (struct qbs*)X2;
	//	if(S1->fn - S2->fn >0)return 1;
	//	if(S1->fn - S2->fn <0)return -1;		
	//	if(S1->fn==S2->fn && S1->fp==S2->fp)
	{
		if (S1->q - S2->q > 0)return 1;
		if (S1->q - S2->q < 0)return -1;
	}
	//	if(S1->fp - S2->fp >0)return 1;
	//	if(S1->fp - S2->fp <0)return -1;	
	return 0;
}
struct gse {
	int cty;
	int gg;
	int n;
};
int compare_cg(const void* X1, const void* X2)
{
	struct gse* S1 = (struct gse*)X1;
	struct gse* S2 = (struct gse*)X2;
	if (S1->gg - S2->gg > 0)return 1;
	if (S1->gg - S2->gg < 0)return -1;
	return 0;
}
struct clust {
	int sta;//start
	int pos;//centr
	int cep;
	double sco;
};

int compare_cur(const void* X1, const void* X2)
{
	// &&&&&&&&&&&&&&&&&&&&&&&&&&&77
	struct clust* S1 = (struct clust*)X1;
	struct clust* S2 = (struct clust*)X2;
	int z1 = S1->sta - S2->sta;
	if (z1 > 0)return 1;
	if (z1 < 0)return -1;
	int x1 = (int)S1->cep;
	int x2 = (int)S2->cep;
	z1 = x1 - x2;
	if (z1 > 0)return 1;
	if (z1 < 0)return -1;
	return 0;
}

void Mix(int* a, int* b)
{
	int buf = *a;
	*a = *b;
	*b = buf;
}
char* TransStr(char* d)
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
// delete symbol 'c' from input string
void DelChar(char* str, char c)
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
int CheckStr(char* d)
{
	int i, len, ret, size;
	ret = size = 0;
	len = strlen(d);
	for (i = 0; i < len; i++)
	{
		if (strchr("atgc", (int)d[i]) != NULL) { continue; }
		else { ret++; }
	}
	return(ret);
}/*
 void GetSost2(char *d, stru5 a)
 {
 int l, i, j;
 int t[5]={10,11,20,21,22};
 for(i=0;i<a.size;i++)sost[i]=0;
 int len=strlen(d), word;
 double x, dx;
 int left0=0, right, left;
 for(i=0;i<LEFT_NUM;i++)left0+=sha[i].part;
 //stru sha[40]={{56,0},{13,1},{89,2},{15,3},{95,4},{6,5},{6,6},{14,7},{10,8},{17,9},{35,10},{44,11},{-1,-1}};
 for(j=0;j<a.size;j++)
 {
 left=left0;
 for(i=0;i<a.wor[j].reg-1;i++)left+=sha[i+LEFT_NUM].part;
 right=left+sha[LEFT_NUM+a.wor[j].reg-1].part;
 word=strlen(a.wor[j].oli);
 //	for(k=LEFT_NUM;sha[k+RIGHT_NUM].num!=-1;k++)
 {
 //	if(a.wor[j].reg-1!=sha[k].num)continue;
 dx=1/(double)word;
 if(LEFT_NUM+a.wor[j].reg!=0)
 {
 x=dx;
 for(i=left-word+1;i<left;i++)
 {
 if(strnicmp(&d[i],a.wor[j].oli,word)==0)
 {
 if(strchr("atgc",a.wor[j].oli[0])!=NULL)
 {
 sost[j]+=x;
 }
 else
 {
 for(l=0;l<5;l++)
 {
 if(i+t[l]<len)
 {
 if(strncmp(&d[i],&d[i+t[l]],word)==0)sost[j]+=x/2;
 }
 if(i-t[l]>=0)
 {
 if(strncmp(&d[i],&d[i-t[l]],word)==0)sost[j]+=x/2;
 }
 }
 }
 }
 }
 x+=dx;
 }
 int left1;
 int right1;
 if(LEFT_NUM+sha[a.wor[j].reg].num!=0)left1=left-word+1;
 else left1=left;
 if(sha[a.wor[j].reg+1].num!=-1)right1=right-word+1;
 else right1=right;
 for(i=left;i<right1;i++)
 {
 if(strnicmp(&d[i],a.wor[j].oli,word)==0)
 {
 if(strchr("atgc",a.wor[j].oli[0])!=NULL)
 {
 sost[j]++;
 }
 else
 {
 for(l=0;l<5;l++)
 {
 if(i+t[l]<len)
 {
 if(strncmp(&d[i],&d[i+t[l]],word)==0)sost[j]+=0.5;
 }
 if(i-t[l]>=0)
 {
 if(strncmp(&d[i],&d[i-t[l]],word)==0)sost[j]+=0.5;
 }
 }
 }
 }
 }
 if(sha[a.wor[j].reg+1].num!=-1)
 {
 x=dx;
 for(i=right-word+1;i<right;i++)
 {
 if(strnicmp(&d[i],a.wor[j].oli,word)==0)
 {
 if(strchr("atgc",a.wor[j].oli[0])!=NULL)
 {
 sost[j]+=x;
 }
 else
 {
 for(l=0;l<5;l++)
 {
 if(i+t[l]<len)
 {
 if(strncmp(&d[i],&d[i+t[l]],word)==0)sost[j]+=x/2;
 }
 if(i-t[l]>=0)
 {
 if(strncmp(&d[i],&d[i-t[l]],word)==0)sost[j]+=x/2;
 }
 }
 }
 }
 x+=dx;
 }
 }
 }
 }
 for(j=0;j<a.size;j++)sost[j]/=sha[a.wor[j].reg].part;
 }
 */
int ConvertSym(int& c)
{
	char four[5] = "atgc";
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

//gauss function
double BackLaplace(double u)
{
	const double pi = 3.141592653589;
	const double step = 0.0001;
	double x, ym, y;
	if (u <= 0) return(-1);
	if (u >= 1) return(-1);
	x = y = 0;
	ym = u * sqrt(2 * pi) / 2;
	while (y < ym)
	{
		x += step / 2;
		y += step * exp(-x * x / 2);
		x += step / 2;
	}
	return(x);
	//return 100;
}
int ComplStr(char* d)
{
	char* d1;
	int i, len;
	len = strlen(d);
	if ((d1 = new char[len + 1]) == NULL) return 0;
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
void ReplaceChar(char* str, char c1, char c2)
{
	int i, len = strlen(str);
	for (i = 0; i < len; i++)
	{
		if (str[i] == c1) str[i] = c2;
	}
}

int ReadSeq(char* file, char* file1, int& n, int& len1)
{
	char head[5000];
	int fl = 0, len;
	char symbol;
	int c;
	//	char cyfr[]="0123456789";
	FILE* in, * out;
	len1 = len = n = 0;
	if ((in = fopen(file, "rt")) == NULL)
	{
		printf("Input file %s can't be opened!\n", file);
		return -1;
	}
	symbol = fgetc(in);
	rewind(in);
	if ((out = fopen(file1, "wt")) == NULL)
	{
		printf("Input file %s can't be opened!\n", file1);
		return -1;
	}
	while ((c = fgetc(in)) != -1)
	{
		if ((char)c == symbol)
		{
			if (n != 0)
			{
				fputc('\n', out);
				printf("\b\b\b\b\b\b\b\b\b%9d", len);
			}
			fputc(c, out);
			if (fgets(head, sizeof(head), in) == NULL)exit(1);
			fprintf(out, "%s", head);
			if (len > len1)len1 = len;
			len = 0;
			n++;
			//printf("\nSeq %5d                    ",n);
			if (n % 10 == 0)printf("\b\b\b\b\b\b\b\b\b%9d", n);
			continue;
		}
		if (strchr("\t\n ", c) != NULL)continue;
		if (strchr("ATGCNatgcn", c) != NULL)
		{
			len++;
			if (len % 100000 == 0)printf("\b\b\b\b\b\b\b\b\b%9d", len);
			if (c < 97) c += 32;
			fputc(c, out);
			/*	if(len==3314)
			{
			//getch();
			}*/
			continue;
		}
		else
		{
			if (ConvertSym(c) == -1)
			{
				printf("Unusual base: sequence N %d, letter position %d\n symbol %c\n%s", n, len, c, head);
			}
			//		printf("Ambiguos base %c replaced at random...", (char)c);
			//	fputc(c,out);
			char cn = 'n';
			fputc(cn, out);
		}
	}
	if (len > len1)len1 = len;
	fclose(in);
	fclose(out);
	return 0;
	//printf("\b\b\b\b\b\b\b\b\b%9d", len);
	//printf("\n");
}
int IdeLet(char c, char* alfabet)
{
	int i, ret = -1;
	for (i = 0; i < 4; i++)
	{
		if (c == alfabet[i]) { ret = i; break; }
	}
	return(ret);
}
void GetWords(int word, int size0, int size, char* w0)
{
	int i, j, zna, che, otn;
	for (i = size0; i < size0 + size; i++)
	{
		che = i - size0;
		zna = size / 4;
		memset(s[i].oli, 0, sizeof(s[i].oli));
		for (j = 0; j < word; j++)
		{
			otn = che / zna;
			s[i].oli[j] = w0[otn];
			che -= otn * zna;
			zna /= 4;
		}
		s[i].oli[word] = '\0';
		//printf("%d\t%s\n",i,s[i].oli);
	}
}

int main(int argc, char* argv[])
{
	int ret = 0, len1, i, j, k, m, nc;
	char mess[300], path_fasta[500], filesta[10], fileend[10], genome[10];
	char sitename[120], file1[500], file_thr_fpr[500], file_out_base[500];	
	double thr;
	FILE* out;
	if (argc != 7)
	{
		printf("andy1_mat_long 1file.path 2char genome (hg38 mm10 dm6 at10) 3sitename 4file.thr_fpr 5pval_crit 6file.output_profile ");
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
	char name_chr_zf[25][3] = { "1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25" };
	int n_chr_zf = 25;
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
	int sizelo_zm73[10] = { 301476924,237917468,232245527,242062272,217959525,169407836,176826311,175377492,157038028,149632204 };
	char name_chr_zm[10][3] = { "1","2","3","4","5","6","7","8","9","10" };
	int n_chr_zm = 10;

	strcpy(path_fasta, argv[1]);
	strcpy(genome, argv[2]);
	strcpy(sitename, argv[3]);
	strcpy(file_thr_fpr, argv[4]);
	double pval_crit = atof(argv[5]);
	pval_crit = -log10(pval_crit);
	strcpy(file_out_base, argv[6]);

	int sizelo1[NCHR];//lengths in bp
	int n_chr = 0;
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
		printf("Genome %s is not recognized\n", genome);
		exit(1);
	}
	city sta;
	sta.get_file(sitename);
	{

		FILE* in_thr;
		if ((in_thr = fopen(file_thr_fpr, "rt")) == NULL)
		{
			printf("Output file %s can't be opened!", file_thr_fpr);
			exit(1);
		}
		char dt[200], sfp[50];
		double thr_prev = 2;
		while (fgets(dt, sizeof(dt), in_thr) != NULL)
		{
			int retu = UnderStol(dt, 1, sfp, sizeof(sfp), '\t');
			if (retu == -1)
			{
				printf("Table Thr..FPR Line Error %s", dt);
				exit(1);
			}
			double fpr_cur = atof(sfp);
			if (fpr_cur < pval_crit)
			{
				thr = thr_prev;
				break;
			}
			thr_prev = atof(dt);
		}
		if (thr_prev == 2)
		{
			printf("FPR table %s Error, Program parmameter %f are bad\n", file_thr_fpr, pval_crit);
			exit(1);
		}
		fclose(in_thr);
	}
	char alfabet[5];	
	strcpy(alfabet, "ACGT");
	alfabet[4] = '\0';	
	double p[SEQLEN];
	char dp[2][SEQLEN + 5];	
	int all_pos = 0;
	int rec_pos = 0;
	int olen = sta.len;

	if ((out = fopen(file_out_base, "wt")) == NULL)
	{
		//printf("Input file can't be opened!\n");
		return -4;
	}
	for(nc=0;nc<n_chr;nc++)
	{
		int all_pos1 = 0;
		int rec_pos1 = 0;
		memset(mess, 0, sizeof(mess));
		strcpy(file1, path_fasta);
		strcat(file1, filesta);
		strcat(file1, name_chr[nc]);
		strcat(file1, fileend);
		len1 = sizelo1[nc];
		int n_time;
		if (len1 > SEQLEN)n_time = 1 + len1 / SEQLEN;
		else n_time = 1;
		FILE* in;
		if ((in = fopen(file1, "rt")) == NULL)
		{
			//printf("Input file %s can't be opened!\n", file1);
			return -3;
		}
		long shift = 0;
		fgets(dp[0], (int)SEQLEN, in);
		fprintf(out, ">Seq %s\tSite %s\tThr %f %f\n", genome, sitename, pval_crit, thr);
		//fprintf (out,"%s",dp);
		shift = ftell(in);
		int shift0 = shift;				
		for (int time = 0; time < n_time; time++)
		{
			fseek(in, shift, SEEK_SET);
			fgets(dp[0], (int)SEQLEN + 1, in);
			int len = strlen(dp[0]);
			strcpy(dp[1], dp[0]);
			ComplStr(dp[1]);
			for (i = 0; i < 2; i++)dp[i][len] = '\0';
			int len2 = len - sta.len;
			for (k = 0; k <= len2; k++)
			{
				char d2[50];
				memset(d2, '\0', sizeof(d2));
				p[k] = -1000;
				double sco2 = 0, sco[2] = { 0,0 };
				for (m = 0; m < 2; m++)
				{
					if (m == 0)strncpy(d2, &dp[0][k], olen);
					else strncpy(d2, &dp[1][len2 - k], olen);
					d2[olen] = '\0';
					if (strchr(d2, 'n') != 0)
					{
						sco2 = -1;
						break;
					}
					for (j = 0; j < sta.size; j++)
					{
						int rlenj = (sta.tot[j].end - sta.tot[j].sta + 1);
						double fm = 0;
						for (i = sta.tot[j].sta; i <= sta.tot[j].end; i++)
						{
							int cod = 4 * IdeLet(d2[i], alfabet) + IdeLet(d2[i + 1], alfabet);
							if (sta.tot[j].num == cod) { fm++; }
						}
						if (fm != 0)
						{
							fm /= rlenj;
							sco[m] += sta.tot[j].buf * fm;
						}
					}
				}
				if (sco2 == 0)
				{
					all_pos1++;
					sco2 = Max(sco[0], sco[1]);
					p[k] = (sco2 - sta.min) / sta.raz;
					if (p[k] >= thr)
					{
						char ori;
						if (sco[0] >= sco[1])ori = '+';
						else ori = '-';
						fprintf(out, "%d\t%.18f\t%c\t%s\n", k + 1, p[k], ori, d2);
						rec_pos1++;
					}					
				}
			}
			if ((time + 1) % 100 == 0)printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%10d\t%d", (time + 1) * SEQLEN, rec_pos1);
			shift += (SEQLEN - sta.len + 1);
		}		
		fclose(in);
		printf("Chr%s %d positions %d sites\n", name_chr[nc], all_pos1, rec_pos1);
		all_pos += all_pos1;
		rec_pos += rec_pos1;
	}
	fclose(out);
	char recfile[500];
	strcpy(recfile, file_out_base);
	strcat(recfile, "_rec_pos.txt");
	if ((out = fopen(recfile, "wt")) == NULL)
	{
		printf("Input file can't be opened!\n");
		exit(1);
	}
	fprintf(out, "%s_(%s,Th %.3f)\t%f\t%d\t%d\t\n", genome, sitename, thr, (double)rec_pos / all_pos, rec_pos, all_pos);
	//	printf("%s_(%s,%d,Z %.3f Th %.3f)\t%.2f\t%d\t%d\t%.2f\t%d\t%d\t\n",argv[1],sitename,cmpl,conf_level,1-conf_level*sigma,(double)rec_seq/nseq,rec_seq,nseq,(double)rec_pos/all_pos,rec_pos,all_pos);	
	return 0;
}

