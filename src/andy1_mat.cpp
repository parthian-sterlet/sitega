#define _CRT_SECURE_NO_WARNINGS

#include  <stdio.h>
#include  <stdlib.h>
#include  <string.h>
#include  <ctype.h>
#include  <math.h>
#include  <time.h>
#define DIM 80

int StrNStr(char *str, char c, int n)
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
void DelHole(char *str)
{
	char *hole;
	hole = strstr(str, "\n");
	if (hole != NULL) *hole = 0;
}
int UnderStol(char *str, int nstol, char *ret, size_t size, char sep)
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
double sost[DIM];
#define Min(a,b) ((a)>(b))? (b):(a);
#define Max(a,b) ((a)>(b))? (a):(b);
double sost1[DIM];
double uw[DIM][DIM];
char **seq;
char **hseq;
int legr[DIM], regr[DIM], tssl[DIM], tssr[DIM];
char symbol;

struct due {
	double buf;
	int sta;
	int end;
	int num;
	void get_copy(due *a);
	void print_all(void);
};
void due::get_copy(due *a)
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
	void get_copy(city *a);
	int get_file(char *file);
	void sort_all(void);
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
void city::get_copy(city *a)
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
//#include "andy1_models.h"
//#include "andy1_get_model.h"

struct qbs {
	double q;
	double fn;
	int fp;
};
int compare_q(const void *X1, const void *X2)
{
	struct qbs *S1 = (struct qbs *)X1;
	struct qbs *S2 = (struct qbs *)X2;
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
int compare_cg(const void *X1, const void *X2)
{
	struct gse *S1 = (struct gse *)X1;
	struct gse *S2 = (struct gse *)X2;
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

int compare_cur(const void *X1, const void *X2)
{
	// &&&&&&&&&&&&&&&&&&&&&&&&&&&77
	struct clust *S1 = (struct clust *)X1;
	struct clust *S2 = (struct clust *)X2;
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

void Mix(int *a, int *b)
{
	int buf = *a;
	*a = *b;
	*b = buf;
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
// delete symbol 'c' from input string
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
int CheckStr(char *d)
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
int ConvertSym(int &c)
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
int ComplStr(char *d)
{
	char *d1;
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
void ReplaceChar(char *str, char c1, char c2)
{
	int i, len = strlen(str);
	for (i = 0; i < len; i++)
	{
		if (str[i] == c1) str[i] = c2;
	}
}

void ReadSeq(char *file, char *file1, int &n, int &len1)
{
	char head[5000];
	int fl = 0, len;
	char symbol;
	int c;	
	FILE  *in, *out;
	len1 = len = n = 0;
	if ((in = fopen(file, "rt")) == NULL)
	{
		printf("Input file %s can't be opened!\n", file);
		exit(1);
	}
	symbol = fgetc(in);
	rewind(in);
	if ((out = fopen(file1, "wt")) == NULL)
	{
		printf("Input file %s can't be opened!\n", file1);
		exit(1);
	}
	while ((c = fgetc(in)) != -1)
	{
		if ((char)c == symbol)
		{
			if (n != 0)
			{
				fputc('\n', out);
				//printf("\b\b\b\b\b\b\b\b\b%9d", len);
			}
			fputc(c, out);
			if (fgets(head, sizeof(head), in) == NULL)exit(1);
			fprintf(out, "%s", head);
			if (len > len1)len1 = len;
			len = 0;
			n++;
			//printf("\nSeq %5d                    ",n);
			//if (n % 10 == 0)printf("\b\b\b\b\b\b\b\b\b%9d", n);
			continue;
		}
		if (strchr("\t\n ", c) != NULL)continue;
		if (strchr("ATGCNatgcn", c) != NULL)
		{
			len++;
		//	if (len % 10000 == 0)printf("\b\b\b\b\b\b\b\b\b%9d", len);
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
	//printf("\b\b\b\b\b\b\b\b\b%9d", len);
	//printf("\n");
}
void MinusStr1(int size, int j, int j0, double buf, double b[DIM][DIM])
{
	int k;
	b[j][j0] = 0;
	for (k = j0 + 1; k < size; k++)
	{
		if (b[j0][k] == 0)continue;
		b[j][k] -= buf * b[j0][k];
	}
}
void MinusStr2(int size, int j, int j0, double buf, double b[DIM][DIM])
{
	int k;
	b[j][j0] = 0;
	for (k = j - 1; k >= j0; k--)
	{
		if (b[j0][k] == 0)continue;
		b[j][k] -= buf * b[j0][k];
	}
}
void MinusStr(int size, int j, int j0, double buf, double b[DIM][DIM])
{
	int k;
	for (k = 0; k < size; k++)
	{
		if (b[j0][k] == 0)continue;
		b[j][k] -= buf * b[j0][k];
	}
}
int BackMat(int size)
{
	int i, j;
	double buf, b[DIM][DIM];
	if (size == 1)
	{
		uw[0][0] = 1 / uw[0][0];
		return 1;
	}
	for (j = 0; j < size; j++)
	{
		for (i = 0; i < size; i++)
		{
			b[i][j] = 0;
		}
	}
	for (i = 0; i < size; i++)b[i][i] = 1;
	for (j = 0; j < size - 1; j++)
	{
		for (i = j + 1; i < size; i++)
		{
			if (uw[i][j] == 0)continue;
			buf = uw[i][j] / uw[j][j];
			MinusStr1(size, i, j, buf, uw);
			MinusStr(size, i, j, buf, b);
			if (fabs(buf) < 1e-035)
			{
				printf("\n1uwij %g(%d,%d) bij%g(%d,%d)\n", uw[i][j], i, j, b[i][j], i, j);
				printf("%g(%d)\n", uw[j][j], j);
				return -1;
				//
			}
			uw[i][j] = 0;
		}
	}
	for (j = size - 1; j > 0; j--)
	{
		for (i = j - 1; i >= 0; i--)
		{
			if (uw[i][j] == 0)continue;
			buf = uw[i][j] / uw[j][j];
			MinusStr2(size, i, j, buf, uw);
			MinusStr(size, i, j, buf, b);
			if (fabs(buf) < 1e-035)
			{
				printf("\n2uwij %g(%d,%d) bij%g(%d,%d)\n", uw[i][j], i, j, b[i][j], i, j);
				printf("%g(%d)\n", uw[j][j], j);
				//
				return -1;
			}
			uw[i][j] = 0;
		}
	}
	double mul = 1;
	for (i = 0; i < size; i++)
	{
		buf = uw[i][i];
		mul *= buf;
		if (fabs(buf) < 1e-030)
		{
			printf("\n3 %g %d\n", uw[i][i], i);
			return -1;
		}
		for (j = 0; j < size; j++)
		{
			uw[i][j] /= buf;
			b[i][j] /= buf;
		}
	}
	if (fabs(mul) < 1e-30)return -1;
	for (i = 0; i < size; i++)
	{
		for (j = 0; j < size; j++)
		{
			uw[i][j] = b[i][j];
		}
	}
	return 1;
}
int IdeLet(char c, char *alfabet)
{
	int i, ret = -1;
	for (i = 0; i < 4; i++)
	{
		if (c == alfabet[i]) { ret = i; break; }
	}
	return(ret);
}
//extern
double Rmah(char *d1, char *d2, city a, char *alfabet)
{
	int k, n, i, j;
	double ret = 0, fn[DIM], fm[DIM];
	for (k = 0; k < a.size; k++)
	{
		int rlenk = (a.tot[k].end - a.tot[k].sta + 1);
		fn[k] = 0;
		for (n = a.tot[k].sta; n <= a.tot[k].end; n++)
		{
			int cod = 4 * IdeLet(d1[n], alfabet) + IdeLet(d1[n + 1], alfabet);
			if (a.tot[k].num == cod) { fn[k]++; break; }
		}
		fn[k] /= rlenk;
		fm[k] = 0;
		for (n = a.tot[k].sta; n <= a.tot[k].end; n++)
		{
			int cod = 4 * IdeLet(d2[n], alfabet) + IdeLet(d2[n + 1], alfabet);
			if (a.tot[k].num == cod) { fm[k]++; break; }
		}
		fm[k] /= rlenk;
	}
	for (i = 0; i < a.size; i++)
	{
		double dret = 0;
		for (j = 0; j < a.size; j++)
		{
			dret += uw[i][j] * (fn[j] - fm[j]);
		}
		ret += dret * (fn[i] - fm[i]);
	}
	return ret;
}
int Fun(char* d, char* head, char* mess, city* sta, double* p, int nseq, char* alfabet)
{
	int i, j, t, err, ret, len;
	char sitename[300];
	strcpy(sitename, sta->site);	
	len = strlen(d);
	FILE* out;
	char addsite[300];
	strcpy(addsite, sitename);
	char file3[300];
	strcpy(file3, "freq_");
	strcat(file3, addsite);
	sta->sort_all();
	int iend = sta->size - 1;	
	static int m = 0;
	if (m == 0)
	{
		if ((out = fopen(file3, "wt")) == NULL)
		{
			printf("Input file %s can't be opened!\n", file3);
			exit(1);
		}
		for (i = 0; i < sta->size; i++)sost1[i] = 0;
		for (i = 0; i < sta->size; i++)
		{
			for (j = 0; j < sta->size; j++)
			{
				uw[i][j] = 0;
			}
		}
		for (j = 0; j < sta->size; j++)tssl[j] = tssr[j] = 0;
		// int tss0=-1;
		int tss0 = 11;
		if (strstr(sitename, "rfam") != NULL)tss0 = -1;
		int left = 1;
		{
			for (i = 0; i < sta->size; i++)
			{
				legr[i] = left + sta->tot[i].sta;
				regr[i] = left + sta->tot[i].end;
				{
					if (legr[i] > tss0)tssl[i] = -tss0;
					else tssl[i] = -(tss0 + 1);
					if (regr[i] > tss0)tssr[i] = -tss0;
					else tssr[i] = -(tss0 + 1);
				}
				fprintf(out, "[%d;%d] %s", legr[i], regr[i], s[sta->tot[i].num].oli);
				if (i != iend)fprintf(out, "\t");
			}
		}
		fprintf(out, "\n");
		fclose(out);
	}
	if (sta->len > len)
	{
		strcpy(mess, "Sequence too short...");
		return(-1);
	}
	err = CheckStr(d);
	if (err != 0)
	{
		{ strcpy(mess, "All ambiguos bases replaced at random..."); ret = 1; }
	}
	else ret = 1;
	int k, n;		
	{
		if ((out = fopen(file3, "at")) == NULL)
		{
			printf("Input file %s can't be opened!\n", file3);
			exit(1);
		}
		m++;				
		for (k = 0; k < sta->size; k++)
		{
			int rlenk = (sta->tot[k].end - sta->tot[k].sta + 1);
			double fm = 0;
			for (n = sta->tot[k].sta; n <= sta->tot[k].end; n++)
			{
				int cod = 4 * IdeLet(d[n], alfabet) + IdeLet(d[n + 1], alfabet);
				if (sta->tot[k].num == cod)fm++;
			}
			fm /= rlenk;
			for (n = 0; n < sta->size; n++)
			{
				if (k == n)uw[k][n] += fm * fm;
				else
				{
					int rlenn = (sta->tot[n].end - sta->tot[n].sta + 1);
					double fn = 0;
					for (t = sta->tot[n].sta; t <= sta->tot[n].end; t++)
					{
						int cod = 4 * IdeLet(d[t], alfabet) + IdeLet(d[t + 1], alfabet);
						if (sta->tot[n].num == cod)fn++;
					}
					fn /= rlenn;
					uw[k][n] += fm * fn;
				}
			}
			sost1[k] += fm;
			fprintf(out, "%f", fm);
			if (k != iend)fprintf(out, "\t");
		}
		fprintf(out, "\n");
		fclose(out);
	}	
	return(ret);
}
void GetWords(int word, int size0, int size, char *w0)
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
void Mix(char *a, char *b)
{
	char buf = *a;
	*a = *b;
	*b = buf;
}
void BigMix1(char *d)
{
	int r;
	int len = strlen(d);
	for (r = 0; r < len - 1; r++) Mix(&d[r], &d[1 + r + (rand() % (len - 1 - r))]);
}
int main(int argc, char *argv[])
{
	char *d, head[500];	
	int n, ret = 0, len, i, k, j, t;
	char mess[300];
	char sitename[120], file1[500], file_thr_fpr[500], file_test[500], file_out_base[500];
	int rec_pos = 0, all_pos = 0, rec_seq = 0;
	double *p, thr;
	FILE  *out;
	if (argc != 7)
	{
		printf("andy1.exe 1file.seq  2sitename 3file thr_fpr 4site_desc(0 no 1 yes) 5pval_crit 6file output_profile ");//7seq_head 8print_pos 9site_desc 10bit
		exit(1);
	}
	strcpy(file_test, argv[1]);
	strcpy(sitename, argv[2]);
	strcpy(file_thr_fpr, argv[3]);
	double pval_crit = atof(argv[5]);
	strcpy(file_out_base, argv[6]);
	int site_desc = atoi(argv[4]);
	pval_crit = -log10(pval_crit);
	{

		FILE *in_thr;
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
	strcpy(alfabet, "acgt");
	alfabet[4] = '\0';

	GetWords(2, 0, 16, alfabet);
	memset(mess, 0, sizeof(mess));
	strcpy(file1, file_out_base);
	strcat(file1, "_strings1.txt");
	int nseq = 0;	
	FILE *in;	
	ReadSeq(file_test, file1, nseq, len);	
	d = new char[len + 2];
	if (d == NULL) { puts("Out of memory..."); exit(1); }	
	p = new double[len];
	if (p == NULL) { puts("Out of memory..."); exit(1); }
	char** d1;
	d1 = new char* [2];
	if (d1 == NULL)
	{
		puts("Out of memory...");
		return(-1);
	}
	for (t = 0; t < 2; t++)
	{
		if ((d1[t] = new char[len + 2]) == NULL)
		{
			puts("Out of memory...");
			return(-1);
		}
	}
	if ((in = fopen(file1, "rt")) == NULL)
	{
		printf("Input file %s can't be opened!\n", file1);
		//getch();
		exit(1);
	}
	FILE *out_best;
	char file_best_score[500];
	memset(file_best_score, 0, sizeof(file_best_score));
	strcpy(file_best_score, file_test);
	strcat(file_best_score, "_bestscosg");
	if ((out_best = fopen(file_best_score, "wt")) == NULL)
	{
		printf("Input file can't be opened!\n");
		exit(1);
	}
	FILE *out_nsite;
	char file_nsite[500];
	memset(file_nsite, 0, sizeof(file_nsite));
	strcpy(file_nsite, file_test);
	strcat(file_nsite, "_nsit");
	if ((out_nsite = fopen(file_nsite, "wt")) == NULL)
	{
		printf("Input file can't be opened!\n");
		exit(1);
	}
	
	if ((out = fopen(file_out_base, "wt")) == NULL)
	{
		printf("Input file can't be opened!\n");
		exit(1);
	}
	city sta;
	sta.get_file(sitename);
	int olen = sta.len;	
	for (n = 0; n < nseq; n++)
	{		
		if (n % 100 == 0)
		{
			printf("\b\b\b\b\b%5d", n + 1);
		}
		if (fgets(head, sizeof(head), in) == NULL)
		{
			printf("Head String %d reading error!", n);
			//getch();
			exit(1);
		}
		DelChar(head, '\n');
		fprintf(out, "%s\n", head);
		if (fgets(d, len + 2, in) == NULL)
		{
			printf("String %d reading error!", n);
			//getch();
			exit(1);
		}
		DelChar(d, '\n');
		TransStr(d);
		int lens = strlen(d);
		for (i = 0; i < lens; i++)p[i] = -1000;		
		if (site_desc == 1)
		{
			ret = Fun(d, head, mess, &sta, p, nseq, alfabet);
			if (ret != 1)
			{
				printf("Fun ret error %s", mess);
			}
		}
		for (t = 0; t < 2; t++)strcpy(d1[t], d);
		ComplStr(d1[1]);
		for (t = 0; t < 2; t++)d1[t][lens] = '\0';
		int len2 = lens - olen;		
		double sco_best=-1000;
		int nsite = 0;
		for (k = 0; k <= len2; k++)
		{
			char d2[50];
			memset(d2, '\0', sizeof(d2));
			p[k] = -1000;
			double sco2 = 0, sco[2] = { 0,0 };
			for (t = 0; t < 2; t++)
			{
				if (t == 0)strncpy(d2, &d1[0][k],olen);
				else strncpy(d2, &d1[1][len2 - k],olen);
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
						sco[t] += sta.tot[j].buf * fm;
					}
				}
			}
			if (sco2 == 0)
			{
				all_pos++;
				sco2 = Max(sco[0], sco[1]);
				p[k] = (sco2 - sta.min) / sta.raz;
				if (p[k] >= thr)
				{
					char ori;
					if (sco[0] >= sco[1])ori = '+';
					else ori = '-';
					fprintf(out,"%d\t%.18f\t%c\t%s\n", k+1, p[k],ori,d2);
					nsite++;
				}
				if (p[k] >= sco_best)sco_best = p[k];
			}
		}		
		fprintf(out_best, "%f\n", sco_best);
		fprintf(out_nsite, "%d\n", nsite);		
		if (nsite > 0)
		{
			rec_seq++;
			rec_pos += nsite;
		}
	}
	fclose(in);
	fclose(out);
	fclose(out_best);	
	fclose(out_nsite);
	char recfile[500];
	strcpy(recfile, file_out_base);
	strcat(recfile, "_rec_pos.txt");

	if ((out = fopen(recfile, "at")) == NULL)
	{
		printf("Input file can't be opened!\n");
		exit(1);
	}
	fprintf(out, "%s\t%s\t%f\t%f\t%d\t%d\t%f\t%d\t%d\t\n", file_test, sitename, thr, (double)rec_seq / nseq, rec_seq, nseq, (double)rec_pos / all_pos, rec_pos, all_pos);
	printf("%s\t%s\t%f\t%f\t%d\t%d\t%f\t%d\t%d\t\n", file_test, sitename, thr, (double)rec_seq / nseq, rec_seq, nseq, (double)rec_pos / all_pos, rec_pos, all_pos);
	fclose(out);
	delete[] d;
	for (t = 0; t < 2; t++)delete[] d1[t];
	delete[] d1;		
	delete[] p;
	return ret;
}

