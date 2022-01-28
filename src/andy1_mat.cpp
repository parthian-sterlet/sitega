#define _CRT_SECURE_NO_WARNINGS

#include  <stdio.h>
#include  <stdlib.h>
#include  <string.h>
#include  <ctype.h>
#include  <math.h>
#include  <time.h>
#define DIM 205
#define NTRAIN 2
#define NCORCOEF 1000

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
//#include "andy1.h"
#define TEN 10000
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
	double c;
	double std;
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
	c = atof(d);
	std = 0.05;
	char sep = '\t', s[20];
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
	//	char cyfr[]="0123456789";
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
			if (len % 10000 == 0)printf("\b\b\b\b\b\b\b\b\b%9d", len);
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
int Fun(char *d, char *head, char *mess, char *file, double *p, int site_desc, int &len0, int end, int nseq, int cmpl, char *alfabet)//double score was 5th argument
{
	int i, j, err, ret, len;
	city sta;

	sta.get_file(file);
	/*if (GetFun(sitename, &sta) == -1)
	{
		printf("Site %s function not found!", sitename);
		exit(1);
	}*/
	char sitename[300];
	strcpy(sitename, sta.site);
	len0 = sta.len;
	len = strlen(d);
	FILE *out;
	char tatabuf[10], addsite[300];
	strcpy(addsite, sitename);
//	sprintf(tatabuf, "%d", sta.len);
//	strcat(addsite, "_");
//	strcat(addsite, tatabuf);
	char file3[300];
	strcpy(file3, "freq_");
	strcat(file3, addsite);
	if (site_desc == 1)sta.sort_all();
	int iend = sta.size - 1;
	static int m = 0;
	if (site_desc == 1 && m == 0)
	{
		//sta.sort_all();		
		if ((out = fopen(file3, "wt")) == NULL)
		{
			printf("Input file %s can't be opened!\n", file3);
			exit(1);
		}
	//	fprintf(out, "\t\t");
		for (i = 0; i < sta.size; i++)sost1[i] = 0;
		for (i = 0; i < sta.size; i++)
		{
			for (j = 0; j < sta.size; j++)
			{
				uw[i][j] = 0;
			}
		}
		for (j = 0; j < sta.size; j++)tssl[j] = tssr[j] = 0;
		// int tss0=-1;
		int tss0 = 11;
		if (strstr(sitename, "rfam") != NULL)tss0 = -1;
		int left = 1;
		{			
			for (i = 0; i < sta.size; i++)
			{
				legr[i] = left + sta.tot[i].sta;
				regr[i] = left + sta.tot[i].end;
				{
					if (legr[i] > tss0)tssl[i] = -tss0;
					else tssl[i] = -(tss0 + 1);
					if (regr[i] > tss0)tssr[i] = -tss0;
					else tssr[i] = -(tss0 + 1);
				}
				fprintf(out, "[%d;%d] %s", legr[i], regr[i], s[sta.tot[i].num].oli);
				if(i!=iend)fprintf(out, "\t");
			}
		}
		fprintf(out, "\n");		
		fclose(out);
	}
	if (sta.len > len)
	{
		strcpy(mess, "Sequence too short...");
		return(-1);
	}
	err = CheckStr(d);
	if (err != 0)
	{
		{strcpy(mess, "All ambiguos bases replaced at random..."); ret = 1; }
	}
	else ret = 1;
	int k, n, t;
	char *d1;
	if ((d1 = new char[len + 1]) == NULL)
	{
		strcpy(mess, "Not enough memory....");
		return(-1);
	}
	for (k = 0; k <= len - sta.len; k++)
	{
		if ((k + 1) % 100000 == 0)
		{
			printf("\b\b\b\b\b\b\b\b\b%9d", k);
		}
		p[k] = 0;
		{
			memset(d1, 0, len+1);
			for (j = 0; j < sta.len; j++)d1[j] = d[k + j]; d1[sta.len] = '\0';
			{
				if (strchr(d1, 'n') != 0)
				{
					p[k] = -1000;
					continue;
				}
			}
			double sco = sta.c;
			for (j = 0; j < sta.size; j++)
			{
				int rlenj = (sta.tot[j].end - sta.tot[j].sta + 1);
				double fm = 0;
				for (i = sta.tot[j].sta; i <= sta.tot[j].end; i++)
				{
					int cod = 4 * IdeLet(d1[i], alfabet) + IdeLet(d1[i + 1], alfabet);
					if (sta.tot[j].num == cod) { fm++; }
				}
				if (fm != 0)
				{
					fm /= rlenj;
					sco += sta.tot[j].buf*fm;
				}
			}
			p[k] = 1 - fabs(sco - 1);
			if (p[k] < -1)p[k] = -1;
		}
	}
	if (site_desc == 1)
	{
		printf("%d\n", m + 1);
		if ((out = fopen(file3, "at")) == NULL)
		{
			printf("Input file %s can't be opened!\n", file3);
			exit(1);
		}
		m++;
		memset(d1, 0, len + 1);
		strcpy(d1, d);
		for (k = 0; k < sta.size; k++)
		{
			int rlenk = (sta.tot[k].end - sta.tot[k].sta + 1);
			double fm = 0;
			for (n = sta.tot[k].sta; n <= sta.tot[k].end; n++)
			{
				int cod = 4 * IdeLet(d1[n], alfabet) + IdeLet(d1[n + 1], alfabet);
				if (sta.tot[k].num == cod)fm++;
			}
			fm /= rlenk;
			for (n = 0; n < sta.size; n++)
			{
				if (k == n)uw[k][n] += fm * fm;
				else
				{
					int rlenn = (sta.tot[n].end - sta.tot[n].sta + 1);
					double fn = 0;
					for (t = sta.tot[n].sta; t <= sta.tot[n].end; t++)
					{
						int cod = 4 * IdeLet(d1[t], alfabet) + IdeLet(d1[t + 1], alfabet);
						if (sta.tot[n].num == cod)fn++;
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
	delete[] d1;
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
	char *d, head[5000];
	char head1[5000];
	int n, ret = 0, len, i, k, len0;
	char mess[300];
	char sitename[120], file1[500], file_thr_fpr[500], file_test[500], file_out_base[500];
	int cmpl, rec_pos = 0, all_pos = 0, rec_seq = 0;
	double *p, thr;
	FILE  *out, *out1;
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
	int bit = 300;
	int site_desc = atoi(argv[4]);
	int head_pr = 1;//atoi(argv[7]);
	int pos_pr = 1;//atoi(argv[8]);
	cmpl = 2;// atoi(argv[5]);	
	{
		
		FILE *in_thr;
		if ((in_thr = fopen(file_thr_fpr, "rt")) == NULL)
		{
			printf("Output file %s can't be opened!", file_thr_fpr);
			exit(1);
		}
		char dt[200], sfp[50];
		double thr_prev=2;
		while (fgets(dt, sizeof(dt), in_thr) != NULL)
		{			
			int retu = UnderStol(dt, 1, sfp, sizeof(sfp), '\t');
			if (retu == -1)
			{
				printf("Table Thr..FPR Line Error %s", dt);
				exit(1);
			}
			double fpr_cur = atof(sfp);
			if (fpr_cur > pval_crit)
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
	double p_zero = 0.9, dp_zero = 1 - p_zero, step_zero = dp_zero / TEN;
	strcpy(alfabet, "acgt");
	//if(strstr(sitename,"chipseq")!=NULL)strcpy(alfabet,"acgt");
	//else strcpy(alfabet,"atgc");
	alfabet[4] = '\0';

	int cnt_mode = 0;
	
	int cmpl2[2] = { 0,1 };
	int cmpl1;
	if (cmpl == 0)cmpl2[1] = -1;
	if (cmpl == 1)cmpl2[0] = -1;
	GetWords(2, 0, 16, alfabet);
	memset(mess, 0, sizeof(mess));
	strcpy(file1, file_out_base);		
	strcat(file1, "_strings1.txt");
	struct thresh {
		int sit;
		int seq;
	};
	thresh den[TEN];
	int nseq = 0;
	int vmix = 100;
	int print_size = 32;
	FILE *in;
	int len1 = 0, len2 = 0;
	ReadSeq(file_test, file1, nseq, len1);
	len = len1;
	d = new char[len + 2];
	if (d == NULL) { puts("Out of memory..."); exit(1); }
	//memset(d,0,sizeof(d));
	p = new double[len];
	if (p == NULL) { puts("Out of memory..."); exit(1); }	
	int nseq1 = nseq;	
	nseq = nseq1;
	if ((in = fopen(file1, "rt")) == NULL)
	{
		printf("Input file %s can't be opened!\n", file1);
		//getch();
		exit(1);
	}
	char plfile[500];
	char plfile0[500];
	strcpy(plfile, file_out_base);
	strcat(plfile, ".dns");
	if ((out1 = fopen(plfile, "wt")) == NULL)
	{
		printf("Input file can't be opened!\n");
		exit(1);
	}
	fclose(out1);
	strcpy(plfile0, file_out_base);
	strcat(plfile0, ".lst");
	if ((out1 = fopen(plfile0, "wt")) == NULL)
	{
		printf("Input file can't be opened!\n");
		exit(1);
	}
	fclose(out1);
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
	for (i = 0; i < TEN; i++)den[i].sit = den[i].seq = 0;
	int *pl, *alive;
	int *pl_all, *alive_all;
	int bit_count, bit_count_all = 1 + len / bit;
	pl = new int[bit_count_all];
	if (pl == NULL) { puts("Out of memory..."); exit(1); }
	pl_all = new int[bit_count_all];
	if (pl_all == NULL) { puts("Out of memory..."); exit(1); }
	alive = new int[bit_count_all];
	if (alive == NULL) { puts("Out of memory..."); exit(1); }
	alive_all = new int[bit_count_all];
	if (alive_all == NULL) { puts("Out of memory..."); exit(1); }
	for (i = 0; i < bit_count_all; i++)
	{
		pl_all[i] = 0;
		alive_all[i] = 0;
	}
	if ((out = fopen(file_out_base, "wt")) == NULL)
	{
		printf("Input file can't be opened!\n");
		exit(1);
	}
	for (n = 0; n < nseq; n++)
	{
		//printf("%d\n",n);
		//printf("\tSeq %d\n",n+1);
		int rec_pos_seq[2] = { 0,0 };
		for (i = 0; i < bit_count_all; i++)
		{
			pl[i] = 0;
			alive[i] = -1;
		}
		int fl_n = 0;
		if (fgets(head, sizeof(head), in) == NULL)
		{
			printf("Head String %d reading error!", n);
			//getch();
			exit(1);
		}
		DelChar(head, '\n');
		strcpy(head1, head);
		//ReplaceChar(head1,'\t', ' ');
		if (fgets(d, len + 2, in) == NULL)
		{
			printf("String %d reading error!", n);
			//getch();
			exit(1);
		}
		DelChar(d, '\n');
		TransStr(d);
		int lens = strlen(d);
		int bit_count_here = 1 + lens / bit;
		int dall_pos = 0;
		int rec_pos_was = rec_pos;
		clust *cur[2];
		int p_score_seq = -5;
		double score_best = -5;
		int pos_best = -1;
		char cep_best = '+';
		int add_flanks = 100;
		for (cmpl1 = 0; cmpl1 < 2; cmpl1++)
		{
			rec_pos_seq[cmpl1] = 0;
			//	printf("%d\t%d\n",n,cmpl1);
			if (cmpl2[cmpl1] == -1)continue;
			if (cmpl2[cmpl1] == 1) if (ComplStr(d) != 1) { puts("Out of memory..."); exit(1); }
			//len=strlen(d);
			//printf("\n%d               ",len);		
			if (n % 100 == 0)printf("\b\b\b\b\b\b%6d", n);
			for (i = 0; i < lens; i++)p[i] = -1000;
			int share;
			if (nseq == 1)share = 1;
			else share = n / (nseq - 1);
			ret = Fun(d, head, mess, sitename, p, site_desc, len0, share, nseq, cmpl, alfabet);
			if (ret != 1)
			{
				printf("Fun ret error %s", mess);
			}
			int lenp = lens - len0 + 1;
			int half0 = len0 / 2;
			//int dpos=(len0+1)%2;
			int dprint = Max(0, (len0 - print_size) / 2);
			char dir[] = "+-";
			/*		if(head_pr==1)
					{
						fprintf(out,"\t%s\t%c\t%s\n",sitename,dir[cmpl2[cmpl1]],head);
					}*/
			if (ret != -1)
			{
				bit_count = 0;
				int dircep = (1 - 2 * cmpl2[cmpl1]);
				{
					int sst, sen;
					if (cmpl2[cmpl1] == 0)
					{
						sst = 0;
						sen = lenp;
					}
					else
					{
						sen = -1;
						sst = lenp - 1;
					}
					bit_count = half0 / bit;
					//if(empty_pos!=0){double zero=0;for(i=1;i<i_sta;i++)fprintf(out,"%d\t%f\n",i,zero);}								
					i = sst;
					do
					{
						if (p[i] != -1000)
						{
							dall_pos++;
							if (alive[bit_count] == -1)
							{
								alive[bit_count] = pl[bit_count] = 0;
							}
							alive[bit_count]++;
						}
						if (p[i] > score_best)
						{
							score_best = p[i];
							if (cmpl1 == 0)
							{
								cep_best = '+';
								pos_best = 1 + i + half0;//ibest=pos_best-half0-1;
							}
							else
							{
								cep_best = '-';
								pos_best = lens - i - half0;//ibest=len_pro1-pos_best-half0;
							}
						}
						if ((i + half0) % bit == 0)
						{
							bit_count++;
						}
						//int p_score=(int)(TEN*((p[i]/4+0.75)));// from - 3 to 0
						int p_score;
						if (p[i] < p_zero)p_score = 0;
						else p_score = (int)((p[i] - p_zero) / step_zero);// from - 1 to 1
					//	if(p_score<0)p_score=0;
						if (p_score > p_score_seq)p_score_seq = p_score;
						for (k = 0; k < p_score; k++)
						{
							den[k].sit++;
						}
						if (p[i] >= thr)
						{
							rec_pos_seq[cmpl1]++;
							pl[bit_count]++;
							//printf("%d: %d\t",i,pl[bit_count]);//add
							rec_pos++;
							if (fl_n == 0)
							{
								rec_seq++;
								fl_n = 1;
							}
						}
						else// && empty_pos==0
						{
							i += dircep;
							continue;
						}
						int posc;
						if (cmpl2[cmpl1] == 0)posc = 1 + i + half0;
						else posc = lens - i - half0;
						/*
						int posc1,  posc2;
						posc1=i+dprint;
						posc2=i+len0-dprint;
						if(pos_pr==1)fprintf(out,"%d\t",posc);
						fprintf(out,"%f\t",p[i]);
						for(int n1=i;n1<posc1;n1++)fprintf(out,"%c",d[n1]);
						for(n1=posc1;n1<posc2;n1++)fprintf(out,"%c",(char)((int)d[n1]-32));
						for(n1=posc2;n1<i+len0;n1++)fprintf(out,"%c",d[n1]);
						fprintf(out,"\n");	*/
						/*fprintf(outsite, "%s\t%s\t", head1, sitename);
						fprintf(outsite, "%d(%c)\n", posc, dir[cmpl1]);
						int j;
						for (j = -add_flanks; j < 0; j++)
						{
							if (i + j >= 0)fprintf(outsite, "%c", d[i + j]);
							else fprintf(outsite, "%c", 'n');
						}
						for (j = i; j < i + len0; j++)fprintf(outsite, "%c", d[j]);
						int iend = i + len0;
						for (j = 0; j < add_flanks; j++)
						{
							if (iend + j < lens)fprintf(outsite, "%c", d[iend + j]);
							else fprintf(outsite, "%c", 'n');
						}
						fprintf(outsite, "\n");*/
						//int p_score=(int)(TEN*(p[i]-thr)/(1-thr));															
						i += dircep;
					} while (i != sen);
					//if(empty_pos!=0){double zero=0;for(i=len-i_end;i<=len;i++)fprintf(out,"%d\t%f\n",i,zero);}
					//printf("\n");
				}
			}
			int rec_pos_here = rec_pos - rec_pos_was;
			cur[cmpl1] = new clust[rec_pos_here];
			if (cur[cmpl1] == NULL) { puts("Out of memory..."); exit(1); }
			int k1 = 0;
			if (cmpl1 != 0 && rec_pos_seq[0] != 0)
			{
				for (i = 0; i < rec_pos_seq[0]; i++)
				{
					cur[1][k1].cep = cur[0][k1].cep;
					cur[1][k1].sta = cur[0][k1].sta;
					cur[1][k1].pos = cur[0][k1].pos;
					cur[1][k1].sco = cur[0][k1].sco;
					k1++;
				}
			}
			for (i = 0; i < lenp; i++)
			{
				if (p[i] >= thr)
				{
					cur[cmpl1][k1].cep = cmpl2[cmpl1];
					cur[cmpl1][k1].sco = p[i];
					if (cmpl2[cmpl1] == 0)
					{
						cur[cmpl1][k1].pos = 1 + i + half0;
						cur[cmpl1][k1].sta = i;
					}
					else
					{
						cur[cmpl1][k1].pos = lens - i - half0;
						cur[cmpl1][k1].sta = lenp - i - 1;
					}
					k1++;
				}
			}
			if (((cmpl == 2 && cmpl1 == 1) || cmpl != 2) && strstr(file_test, "~") == NULL)// random omitted!
			{
				int j;
				if (cmpl >= 1) if (ComplStr(d) != 1) { puts("Out of memory..."); exit(1); }
				int rec_total = rec_pos_seq[0] + rec_pos_seq[1];
				int cur_in = 1;
				if (cmpl == 0)cur_in = 0;
				char dir[] = "+-";
				//	fprintf(out_nsite,"%d\t%s\t%d\n",n+1,head1,rec_total);	
				fprintf(out_nsite, "%d\n", rec_total);
				//	fprintf(out_best,"%d%c\t%f\t",pos_best,cep_best,score_best);			
				fprintf(out_best, "%.12f", score_best);
				//fprintf(out_best,"\t");			
				/*{
					int n1, posc1, posc2;
					char d3[200];
					int ipos;
					if(cep_best=='+')ipos=pos_best-half0-1;
					else ipos=pos_best+half0-len0;
					posc1=ipos+dprint;
					posc2=ipos+len0-dprint;
					if(posc1<0)posc1=0;
					if(posc2>lens-1)posc2=lens-1;
					memset(d3,0,sizeof(d3));
					for(n1=posc1;n1<posc2;n1++)d3[n1-posc1]=d[n1];
					if(cep_best=='-')if(ComplStr(d3)!=1) {puts("Out of memory...");exit(1);}
					fprintf(out_best,"%s",d3);
				}*/
				fprintf(out_best, "\n");
				if (rec_total > 0)
				{
					FILE *outpos;
					if ((outpos = fopen(plfile, "at")) == NULL)
					{
						printf("Input file can't be opened!\n");
						exit(1);
					}
					int bit0 = half0 / bit;
					if (rec_seq == 1)
					{
						//fprintf(outpos,"%.2f\n#\tHead\tSite\tPos+\tPos-\tSeq+\tSeq-\tSco+\tSco-\t",thr);
						fprintf(outpos, "%.3f\t%s\tSite\tPosition\tStrand\tSequence\tScore\tLength", thr, argv[0]);
						/*	for(i=bit0;i<bit_count_all-bit0-1;i++)
							{
								fprintf(outpos,"\t");
								fprintf(outpos,"%d",bit*i+bit/2);

							}*/
						fprintf(outpos, "\n");
					}
					//fprintf(outpos,"%d\t%s\t%s\t",n+1,head1,sitename);
					for (j = 0; j < rec_total; j++)
					{
						//if(j==rec_pos_seq[0])fprintf(outpos,"\t");
						fprintf(outpos, "%d\t%s\t%s\t", n + 1, head1, sitename);
						fprintf(outpos, "%d\t(%c)\t", cur[cur_in][j].pos, dir[cur[cur_in][j].cep]);
						int n1, posc1, posc2;
						int ipos;
						//int dprint=Max(0,len1/2-5);	
						int dprint = Max(0, len0 / 2 - 5);
						//int dprint=Max(0,len1/2+25);	
						if (cur[cur_in][j].cep == 0)ipos = cur[cur_in][j].pos - half0 - 1;
						else
						{
							//ipos=cur[cur_in][j].pos+half0-len1;
							ipos = cur[cur_in][j].pos + half0 - len0;
						}
						posc1 = cur[cur_in][j].pos - dprint;
						posc2 = cur[cur_in][j].pos + dprint;
						char d3[200];
						memset(d3, 0, sizeof(d3));
						for (n1 = ipos; n1 < posc1; n1++)d3[n1 - ipos] = d[n1];
						for (n1 = posc1; n1 < posc2; n1++)d3[n1 - ipos] = (char)((int)(d[n1] - 32));
						for (n1 = posc2; n1 < ipos + len0; n1++)d3[n1 - ipos] = d[n1];
						//strcpy(d3,dp);
						if (cur[cur_in][j].cep == 1)if (ComplStr(d3) != 1) { puts("Out of memory..."); exit(1); }
						fprintf(outpos, "%s\t", d3);
						fprintf(outpos, "%.6f\t%d\n", cur[cur_in][j].sco, lens);
					}
					fclose(outpos);
					/*
					for(j=0;j<rec_total;j++)
					{
						if(j==rec_pos_seq[0])fprintf(outpos,"\t");
						fprintf(outpos,"%d(%c) ",cur[cur_in][j].pos,dir[cur[cur_in][j].cep]);
					}
					if(rec_pos_seq[1]==0)fprintf(outpos,"\t");
					fprintf(outpos,"\t");
					int dprint=Max(0,len0/2-5);
					for(j=0;j<rec_total;j++)
					{
						if(j==rec_pos_seq[0])fprintf(outpos,"\t");
						int n1, posc1, posc2;
						char d3[200];
						int ipos;
						if(cur[cur_in][j].cep==0)ipos=cur[cur_in][j].pos-half0-1;
						else
						{
							ipos=cur[cur_in][j].pos+half0-len0;
						}
						posc1=ipos+dprint;
						posc2=ipos+len0-dprint;
						if(posc1<0)posc1=0;
						if(posc2>lens-1)posc2=lens-1;
						memset(d3,0,sizeof(d3));
						for(n1=posc1;n1<posc2;n1++)d3[n1-posc1]=d[n1];
						if(cur[cur_in][j].cep==1)if(ComplStr(d3)!=1) {puts("Out of memory...");exit(1);}
						fprintf(outpos,"%s",d3);
						fprintf(outpos," ");
					}
					if(rec_total==rec_pos_seq[0])fprintf(outpos,"\t");
					if(rec_pos_seq[1]==0)fprintf(outpos,"\t");
					for(j=0;j<rec_total;j++)if(cur[cur_in][j].cep==0){fprintf(outpos,"%.3f",cur[cur_in][j].sco);if(j!=rec_total-1)fprintf(outpos,",");}
					fprintf(outpos,"\t");
					for(j=0;j<rec_total;j++)if(cur[cur_in][j].cep==1){fprintf(outpos,"%.3f",cur[cur_in][j].sco);if(j!=rec_total-1)fprintf(outpos,",");}
					fprintf(outpos,"\t");
					*/
					//fprintf(outpos,"\t");				
				}
				{
					FILE *outpos;
					if ((outpos = fopen(plfile0, "at")) == NULL)
					{
						printf("Input file can't be opened!\n");
						exit(1);
					}
					int bit0 = half0 / bit;
					if (n == 0)
					{
						fprintf(outpos, "#\tHead\tSite\tPos+\tPos-\tSco+\tSco-\tSeq+\tSeq-");
						for (i = bit0; i < bit_count_all - bit0 - 1; i++)
						{
							fprintf(outpos, "\t");
							fprintf(outpos, "%d", bit*i + bit / 2);

						}
						fprintf(outpos, "\n");
					}
					fprintf(outpos, "%d\t%s\t%s\t", n + 1, head1, sitename);
					int tss_pos = 3000;
					for (j = 0; j < rec_total; j++)
					{
						if (j == rec_pos_seq[0])fprintf(outpos, "\t");
						fprintf(outpos, "%d(%c) ", tss_pos - cur[cur_in][j].pos, dir[cur[cur_in][j].cep]);
					}
					if (rec_pos_seq[1] == 0)fprintf(outpos, "\t");
					fprintf(outpos, "\t");
					for (j = 0; j < rec_total; j++)
					{
						if (j == rec_pos_seq[0])fprintf(outpos, "\t");
						fprintf(outpos, "%.3f ", cur[cur_in][j].sco);
					}
					if (rec_pos_seq[1] == 0)fprintf(outpos, "\t");
					fprintf(outpos, "\t");
					for (j = 0; j < rec_total; j++)
					{
						if (j == rec_pos_seq[0])fprintf(outpos, "\t");
						int n1, posc1, posc2;
						char d3[200];
						int ipos;
						if (cur[cur_in][j].cep == 0)ipos = cur[cur_in][j].pos - half0 - 1;
						else
						{
							//	ipos=cur[cur_in][j].pos-half0-1-dpos;
							ipos = cur[cur_in][j].pos + half0 - len0;
						}
						posc1 = ipos + dprint;
						posc2 = ipos + len0 - dprint;
						if (posc1 < 0)posc1 = 0;
						if (posc2 > lens - 1)posc2 = lens - 1;
						memset(d3, 0, sizeof(d3));
						for (n1 = posc1; n1 < posc2; n1++)d3[n1 - posc1] = d[n1];
						if (cur[cur_in][j].cep == 1)if (ComplStr(d3) != 1) { puts("Out of memory..."); exit(1); }
						fprintf(outpos, "%s", d3);
						fprintf(outpos, " ");
					}
					if (rec_pos_seq[1] == 0)fprintf(outpos, "\t");
					//fprintf(outpos,"\t");
					fclose(outpos);
				}
				qsort(cur[cmpl1], rec_pos_here, sizeof(cur[cmpl1][0]), compare_cur);
				fprintf(out, "%s\t%s\tSEQ %d\tTHR %f\n", head, argv[2], n + 1, thr);
				{
					for (j = 0; j < rec_total; j++)
					{
						fprintf(out, "%d\t%.18f\t%c\t", cur[cur_in][j].sta, cur[cur_in][j].sco, dir[cur[cur_in][j].cep]);
						int n1, posc1, posc2;
						int dprint = 4;//Max(10,len1/2);
						char d3[200];
						int ipos;
						if (cur[cur_in][j].cep == 0)ipos = cur[cur_in][j].pos - half0 - 1;
						else
						{
							ipos = cur[cur_in][j].pos + half0 - len0;
						}
						//posc1=cur[cur_in][j].pos-dprint;											
						//posc2=cur[cur_in][j].pos+dprint;
						posc1 = cur[cur_in][j].sta;
						posc2 = posc1 + len0;
						if (posc1 < 0)posc1 = 0;
						if (posc2 > lens)posc2 = lens;
						memset(d3, 0, sizeof(d3));
						for (n1 = posc1; n1 < posc2; n1++)d3[n1 - posc1] = d[n1];
						d3[len0] = '\0';
						if (cur[cur_in][j].cep == 1)if (ComplStr(d3) != 1) { puts("Out of memory..."); exit(1); }
						int poscap1 = half0 - dprint, poscap2 = half0 + dprint;
						//if(poscap1<0)poscap1=0;
						//if(poscap2>lens-1)poscap2=lens-1;
						for (n1 = poscap1; n1 < poscap2; n1++)d3[n1] = (char)((int)(d3[n1] - 32));
						fprintf(out, "%s", d3);
						fprintf(out, "\n");
					}
				}
			}
		}
		//	if(i_best!=-1)
		if (p_score_seq != -1)
		{

			for (i = 0; i < p_score_seq; i++)
			{
				den[i].seq++;
			}
		}
		/*if (fl_n != 0)
		{
			FILE *out1;
			if ((out1 = fopen("yes.txt", "at")) == NULL)
			{
				printf("Input file can't be opened!\n");
				exit(1);
			}
			if (cmpl2[cmpl1] == 1)if (ComplStr(d) != 1) { puts("Out of memory..."); exit(1); }
			fprintf(out1, "%s\n%s\n", head1, d);
			fclose(out1);
		}*/
		if (cmpl == 2)dall_pos /= 2;
		all_pos += dall_pos;
		/*	if(bit_count_all>255)
			{
				fprintf(out,"\t%s_(%d,%.2f)\n",file_test,tata,thr);
				for(bit0=0;i<bit_count_here-bit0;i++)
				{
					fprintf(out,"%d",bit*(1+i));
					fprintf(out,"\t");
					if(alive[i]!=-1)fprintf(out,"%d",pl[i]);
					fprintf(out,"\n");
				}
			}
			else*/
		int bit0 = len0 / bit / 2;
		/*
		if(rec_pos_seq[0]+rec_pos_seq[1]!=0)
		{
			if((out1=fopen(plfile,"at"))==NULL)
			{
				printf("Input file can't be opened!\n");
				exit(1);
			}
			//int a1=cmpl2[cmpl1]*(bit_count+1)*bit;
			//int a2=(-2*cmpl2[cmpl1]+1)*bit;
			int alive_all_here=0;
			for(i=bit0;i<bit_count_here-bit0-1;i++)
			{
				fprintf(out1,"\t");
				if(alive[i]!=-1)
				{
					alive_all_here++;
					fprintf(out1,"%d",pl[i]);
				}
			}
			fprintf(out1,"\t%d",alive_all_here*bit);
			fprintf(out1,"\n");
			fclose(out1);
		}*/
		{
			if ((out1 = fopen(plfile0, "at")) == NULL)
			{
				printf("Input file can't be opened!\n");
				exit(1);
			}
			//int a1=cmpl2[cmpl1]*(bit_count+1)*bit;
			//int a2=(-2*cmpl2[cmpl1]+1)*bit;	
			int alive_all_here = 0;
			for (i = bit0; i < bit_count_here - bit0 - 1; i++)
			{
				fprintf(out1, "\t");
				if (alive[i] != -1)
				{
					alive_all_here++;
					fprintf(out1, "%d", pl[i]);
				}
			}
			fprintf(out1, "\t%d", alive_all_here*bit);
			fprintf(out1, "\n");
			fclose(out1);
		}
		/*	if((out=fopen("site_table.txt","at"))==NULL)
			{
				printf("Input file can't be opened!\n");
				exit(1);
			}
			fprintf(out,"%s\t",head1);
			for(i=0;i<lenp;i++)
			{
				fprintf(out,"%f\t",p[i]);
			}
			fprintf(out,"\n");
			fclose (out);
		*/	for (i = bit0; i < bit_count_here - bit0; i++)
		{
			pl_all[i] += pl[i];
			if (cmpl == 2)alive[i] /= 2;
			alive_all[i] += alive[i];
			//printf("%d: %d\t",i,pl_all[i]);//add
		}
		//printf("\n");
	}
	fclose(in);
	fclose(out);
//	fclose(outsite);
	fclose(out_best);
	fclose(out_nsite);
	int bit0 = len0 / bit / 2;
	char dnsfile[500];
	strcpy(dnsfile, file_out_base);
	strcat(dnsfile, "_dns_pos.txt");
	if (cmpl != 1)
	{
		if ((out = fopen(dnsfile, "at")) == NULL)
		{
			printf("Input file can't be opened!\n");
			exit(1);
		}
		fprintf(out, "\t");
		for (i = bit0; i < 1 + bit_count_all - bit0; i++)
		{
			fprintf(out, "%d", bit*(1 + i) - bit / 2);
			fprintf(out, "\t");
		}
		fprintf(out, "\n");
	}
	else
	{
		if ((out = fopen(dnsfile, "at")) == NULL)
		{
			printf("Input file can't be opened!\n");
			exit(1);
		}
	}
	if (cmpl != 2)
	{
		fprintf(out, "SITE(%d)\t", cmpl);
		for (i = bit0; i < 1 + bit_count_all - bit0; i++)
		{
			if (alive_all[i] != 0)fprintf(out, "%d", pl_all[i]);
			fprintf(out, "\t");
		}
		fprintf(out, "\n");
		fprintf(out, "ALIVE(%d)\t", cmpl);
		for (i = bit0; i < 1 + bit_count_all - bit0; i++)
		{
			fprintf(out, "%d", alive_all[i]);
			fprintf(out, "\t");
		}
		fprintf(out, "\n");
	}
	else
	{
		fprintf(out, "%s_DENS(%d)", file_test, cmpl);
		for (i = bit0; i < 1 + bit_count_all - bit0; i++)
		{
			if (alive_all[i] != 0)fprintf(out, "\t%f", (double)pl_all[i] / alive_all[i]);
		}
		fprintf(out, "\n");
	}
	//	delete(p);
	fclose(out);
	char recfile[500];
	strcpy(recfile, file_out_base);
	strcat(recfile, "_rec_pos.txt");

	if ((out = fopen(recfile, "at")) == NULL)
	{
		printf("Input file can't be opened!\n");
		exit(1);
	}
	fprintf(out, "%s_(%s,%d)\t%f\t%f\t%d\t%d\t%f\t%d\t%d\t\n", file_test, sitename, cmpl, thr, (double)rec_seq / nseq, rec_seq, nseq, (double)rec_pos / all_pos, rec_pos, all_pos);
	printf("%s_(%s,%d, Th %f)\t%.2f\t%d\t%d\t%.2f\t%d\t%d\t\n", file_test, sitename, cmpl, thr, (double)rec_seq / nseq, rec_seq, nseq, (double)rec_pos / all_pos, rec_pos, all_pos);
	//	delete(p);
	fclose(out);
	char pltfile[500];
	strcpy(pltfile, file_out_base);
	strcat(pltfile, ".plt");
	if ((out = fopen(pltfile, "at")) == NULL)
	{
		printf("Input file can't be opened!\n");
		exit(1);
	}
	//	fprintf(out,"\t");
		/*
		fprintf(out,"Thr\t");
		fprintf(out,"SEQ_%s_%s\t",file_test,sitename);
		fprintf(out,"SITE_%s_%s\n",file_test,sitename);
		for(i=0;i<TEN;i++)
		{
			fprintf(out,"%.6f\t",thr+(1-thr)*i/TEN);
			fprintf(out,"%d\t",den[i].seq);
			fprintf(out,"%d\n",den[i].sit);
		}*/
	fprintf(out, "SEQ_%s", file_test);
	for (i = TEN - 1; i >= 0; i--)
	{
		fprintf(out, "\t%d", den[i].seq);
	}
	fprintf(out, "\n");
	//	if(strstr(file_test,"~")!=NULL)
	{
		int ten1 = TEN - 1;
		fprintf(out, "SITE_%s", file_test);
		for (i = ten1; i >= 0; i--)
		{
			fprintf(out, "\t%d", den[i].sit);
		}
		fprintf(out, "\n");
		fprintf(out, "Thr_%s_%s", file_test, sitename);
		double thr_max = 1;
		for (i = TEN; i >= 1; i--)
		{
			//fprintf(out,"\t%f",thr+thr_m1*i/ten1);	
			fprintf(out, "\t%.7f", thr_max);
			thr_max -= step_zero;
			//fprintf(out,"\t%f",-3+4*(double)(i)/TEN);	
		}
		fprintf(out, "\n");
	}
	fclose(out);
	delete[] d;
	delete[] p;
	return ret;
}

