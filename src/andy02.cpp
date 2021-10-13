#define _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS

#include  <stdio.h>
#include  <stdlib.h>
#include  <string.h>
//#include  <conio.h>
#include  <math.h>
#include  <time.h>
#include <ctype.h>

#define Min(a,b) ((a)>(b))? (b):(a);
#define Max(a,b) ((a)>(b))? (a):(b);
#define SEQLEN 12000
#define MOTLEN 12 //max LPD length
#define MEGE 40//population size 1st stage
#define ELIT 15//population size 2nd stage
#define NMUT 3
#define NREC 5
#define POPSIZE 200
#define CENT 100

double  uw[POPSIZE][POPSIZE], uw0[POPSIZE];
double df[POPSIZE];
struct ss{
	int num;
	char oli[3];
}s[16];

int compare_qq(const void *X1, const void *X2)
{
	double X = (*(double*)X1 - *(double*)X2);
	if (X>0)return 1;
	if (X<0)return -1;
	return 0;
}
int compare_qq2(const void *X1, const void *X2)
{
	double X = (*(double*)X1 - *(double*)X2);
	if (X>0)return -1;
	if (X<0)return 1;
	return 0;
}
struct uno{
	int sta;
	int end;
	int num;
	void get_copy(uno *a);
	void print_all(void);
};
void uno::get_copy(uno *a)
{
	a->num = num;
	a->sta = sta;
	a->end = end;
};
void uno::print_all(void)
{
	printf("[%d;%d]%s\t", sta, end, s[num].oli);
}
struct town{
	uno tot[POPSIZE];
	int deg[16];
	int size;
	int *pos;// pozicii na4al okon
	int *ori;// DNA strand 0,1
	double mah;
	double std;
	double ave;
	double fit;
	int odg[MOTLEN+2];
	void get_copy(town *a, int nseq, int reg_max);
	void init_rand(int nseq, int *len, int oln, int rsize, int reg_max);
	void init_rand_part(int nseq, int *len, int oln, int nind);
	int init_add(uno last);
	int init_add_split(void);
	void init_zero(int olen);
	void Mix(int *a, int *b);
	void BigMix(int *d, int ln);
	void swap(int n1, int n2);
	int order(int n);
	void print_all(int reg_max, int nseq);
	void print_sta(int reg_max);
	void sort_all(void);
	int sum(int j);
	void fprint_all(char *file, char *add);
	void fprint_allfi(char *file, char *add, int len, double c0, double *buf);
	void fprint_allfi_mat(char *file, char *add, char *name, int len, double c0, double *buf);
	void fprint_seq(char *file, int len, int nseq, char ***seq, double *best_sco);
	int check(int min, int max);
	int mem_in(int nseq);
	void mem_out(void);
} pop[MEGE], det1, det2[2];
int compare_tot(const void *X1, const void *X2)
{
	struct uno *S1 = (struct uno *)X1;
	struct uno *S2 = (struct uno *)X2;
	if (S1->num - S2->num >0)return 1;
	if (S1->num - S2->num <0)return -1;
	if (S1->sta - S2->sta >0)return 1;
	if (S1->sta - S2->sta <0)return -1;
	return 0;
}
int compare_pop(const void *X1, const void *X2)
{
	struct town *S1 = (struct town *)X1;
	struct town *S2 = (struct town *)X2;
	if (S1->fit - S2->fit >0)return -1;
	if (S1->fit - S2->fit <0)return 1;
	return 0;
}
void town::init_zero(int reg_max)
{
	int k;
	for (k = 0; k<16; k++)deg[k] = 0;
	for (k = 0; k<reg_max; k++)odg[k] = 0;
	odg[reg_max] = -1;
	size = 0;
	fit = mah = std = ave = 0;
}
int town::init_add(uno last)
{
	int k;
	int dtype = last.num;
	int sta = last.sta;
	int end = last.end;
	for (k = 0; k<size; k++)
	{
		if (tot[k].num == dtype)
		{
			int sta1 = tot[k].sta;
			int end1 = tot[k].end;
			if (sta >= sta1 && sta <= end1)return -1;
			if (end >= sta1 && end <= end1)return -1;
			if (sta1 >= sta && sta1 <= end)return -1;
			if (end1 >= sta && end1 <= end)return -1;
		}
	}
	tot[size].sta = sta;
	tot[size].end = end;
	tot[size].num = dtype;
	deg[dtype]++;
	odg[end - sta]++;
	size++;
	return 1;
}
int town::init_add_split(void)
{
	int k, r, r1, x, end1;
	do
	{
		r = rand() % size;
	} while (tot[r].end - tot[r].sta == 0);
	for (k = size - 1; k > r; k--)
	{
		tot[k].get_copy(&tot[k + 1]);
	}
	x = tot[r].end - tot[r].sta;
	end1 = rand() % x;
	r1 = r + 1;
	tot[r1].sta = tot[r].sta + end1 + 1;
	tot[r1].end = tot[r].end;
	tot[r1].num = tot[r].num;
	tot[r].end = tot[r].sta + end1;
	deg[tot[r].num]++;
	odg[x]--;
	odg[tot[r1].end - tot[r1].sta]++;
	odg[tot[r].end - tot[r].sta]++;
	size++;
	return 1;
}
void town::sort_all(void)
{
	qsort((void*)tot, size, sizeof(tot[0]), compare_tot);
}
int town::check(int min, int max)
{
	int i, j;
	int ods[MOTLEN];
	int max1 = max - 1;
	for (i = 0; i<max1; i++)ods[i] = 0;
	for (i = 0; i<size; i++)
	{
		ods[tot[i].end - tot[i].sta]++;
	}
	for (i = 0; i<max1; i++)
	{
		if (ods[i] != odg[i])
		{
			printf("odg(%d) %d actual %d", i, odg[i], ods[i]);
			return -1;
		}
	}
	int sum = 0;
	for (i = 0; i<max; i++)
	{
		if (odg[i]<0){ printf("odg %d = %d", i, odg[i]); return -1; }
		sum += odg[i];
	}
	if (sum != size)
	{
		printf("wrong sum odg = %d, size = %d", sum, size);
		return -1;
	}
	for (i = 0; i<size; i++)
	{
		int len = tot[i].end - tot[i].sta;
		if (len<min || len >= max)
		{
			tot[i].print_all();
			printf("\n");
			return -1;
		}
	}
	for (i = 0; i<size; i++)
	{
		for (j = i + 1; j<size; j++)
		{
			if (tot[i].num == tot[j].num)
			{
				if (tot[i].end <= tot[j].sta)continue;
				if (tot[j].end <= tot[i].sta)continue;
				int cat = 0;
				if (tot[i].sta >= tot[j].sta && tot[j].sta <= tot[i].end)cat = 1;
				if (tot[j].sta >= tot[i].sta && tot[i].sta <= tot[j].end)cat = 1;
				if (tot[i].sta >= tot[j].end && tot[j].end <= tot[i].end)cat = 1;
				if (tot[j].sta >= tot[i].end && tot[i].end <= tot[j].end)cat = 1;
				if (cat == 1)
				{
					printf("1st %d\t", i + 1);
					tot[j].print_all();
					printf("\n");
					printf("2nd %d\t", j + 1);
					tot[i].print_all();
					printf("\n");
					return -1;
				}
			}
		}
	}
	return 1;
}
struct town_ext{
	double c0;
	double buf[POPSIZE];
	void get_copy(double c1, double *b1, int size);
} pop_ext;
void town_ext::get_copy(double c1, double *b1, int size)
{
	c0 = c1;
	int i;
	for (i = 0; i<size; i++)buf[i] = b1[i];
}
void MixI(int *a, int *b)
{
	int buf = *a;
	*a = *b;
	*b = buf;
}
void BigMixI(int *d1, int len) // pereme6ivanie stroki
{
	int r;
	for (r = 0; r<len - 1; r++)
	{
		MixI(&d1[r], &d1[1 + r + (rand() % (len - 1 - r))]);
	}
}
void MixC(char *a, char *b)
{
	char buf = *a;
	*a = *b;
	*b = buf;
}
void BigMixC(char *d1, int len) // pereme6ivanie stroki
{
	int r;
	for (r = 0; r<len - 1; r++)
	{
		MixC(&d1[r], &d1[1 + r + (rand() % (len - 1 - r))]);
	}
}
void town::print_all(int reg_max, int nseq)
{
	int i;
	char strand[] = "+-";
	printf("M %f A %f S %f F %f\t", mah, ave, std, fit);
	for (i = 0; i<16; i++)printf("%d ", deg[i]); printf("\tLEN ");
	for (i = 0; i<reg_max; i++)printf("%d ", odg[i]); printf("\t");
	int size1 = size - 1;
	printf(" ");
	for (i = 0; i<size; i++)
	{
		tot[i].print_all();
		if (i == size1)printf("\n");
		// else if(tot[i + 1].num != tot[i].num)printf("\n");
		//	 printf(" ");
	}
	// printf("\n");
	/*
	for(i=0;i<nseq;i++)
	{
	printf("%2d",pos[i]);
	printf("%c",strand[ori[i]]);
	printf(" ");
	}
	printf("\n");*/
}
void town::print_sta(int reg_max)
{
	int i;
	printf("M %f A %f S %f F %f\t", mah, ave, std, fit);
	for (i = 0; i<16; i++)printf("%2d", deg[i]); printf("\t");
	for (i = 0; i<reg_max; i++)printf("%2d ", odg[i]); printf("\n");
}
void town::get_copy(town *a, int nseq, int reg_max)
{
	int i;
	a->size = size;
	a->fit = fit;
	a->ave = ave;
	a->mah = mah;
	a->std = std;
	for (i = 0; i < size; i++)tot[i].get_copy(&a->tot[i]);
	for (i = 0; i < 16; i++)a->deg[i] = deg[i];
	i = 0;
	for (i = 0; i<reg_max; i++)
	{
		a->odg[i] = odg[i];
	}
	for (i = 0; i < nseq; i++)
	{
		a->pos[i] = pos[i];
		a->ori[i] = ori[i];
	}
}
void town::swap(int n1, int n2)
{
	uno buf;
	tot[n1].get_copy(&buf);
	tot[n2].get_copy(&tot[n1]);
	buf.get_copy(&tot[n2]);
}
int town::order(int n)
{
	/*printf("\n");
	printf("old %d\t",n);
	tot[n].print_all();printf("\n");
	print_all();
	*/
	int ret = n;
	int i;
	int di = 0;
	if (n != 0 && n != size - 1)
	{
		if (tot[n].num<tot[n - 1].num)di = -1;
		else if (tot[n].num>tot[n + 1].num)di = 1;
		else
		{
			if (tot[n].num == tot[n - 1].num && tot[n].sta<tot[n - 1].sta)di = -1;
			if (tot[n].num == tot[n + 1].num && tot[n].sta>tot[n + 1].sta)di = 1;
		}
	}
	{
		if (n == 0)
		{
			if (tot[n].num>tot[n + 1].num)di = 1;
			else if (tot[n].num == tot[n + 1].num && tot[n].sta>tot[n + 1].sta)di = 1;
		}
		if (n == size - 1)
		{
			if (tot[n].num<tot[n - 1].num)di = -1;
			else if (tot[n].num == tot[n - 1].num && tot[n].sta<tot[n - 1].sta)di = -1;
		}
	}
	//printf("DIRECT %d\t",di);
	//if(di==0)return -1;
	if (di == -1)
	{
		for (i = 0; i <= n - 1; i++)
		{
			if (tot[n - i - 1].num>tot[n - i].num)
			{
				//			tot[n-i-1].print_all();
				swap(n - i - 1, n - i);
				ret--;
			}
			else if (tot[n - i - 1].num == tot[n - i].num && tot[n - i - 1].sta>tot[n - i].sta)
			{
				//			tot[n-i-1].print_all();
				swap(n - i - 1, n - i);
				ret--;
			}
			else break;
		}
	}
	if (di == 1)
	{
		for (i = 0; i <= size - n - 2; i++)
		{
			if (tot[n + i + 1].num<tot[n + i].num)
			{
				//			tot[n+i+1].print_all();
				swap(n + i + 1, n + i);
				ret++;
			}
			else if (tot[n + i + 1].num == tot[n + i].num && tot[n + i + 1].sta<tot[n + i].sta)
			{
				//			tot[n+i+1].print_all();
				swap(n + i + 1, n + i);
				ret++;
			}
			else break;
		}
	}
	/*printf("\n");
	printf("new %d\t",ret);
	tot[ret].print_all();printf("\n");
	print_all();
	printf("\n");*/
	return ret;
}
void town::init_rand_part(int nseq, int *len, int oln, int nind)
{
	int i;

	fit = 0;
	int oln1 = oln - 1;
	for (i = 0; i < nseq; i++)
	{
		int r = rand() % nind;
		if (r == 0)
		{
			int lenp = len[i] - oln1;
			pos[i] = rand() % lenp;
			ori[i] = rand() % 2;
		}
	}
}
void town::init_rand(int nseq, int *len, int oln, int rsize, int reg_max)
{
	int i, j;

	fit = 0;
	size = rsize;
	int oln1 = oln - 1;
	int oln2 = oln - 2;
	for (j = 0; j < reg_max; j++)odg[j] = 0;
	odg[reg_max] = -1;
	//for (i = 0; i < nseq; i++)printf("%d ",len[i]);
	for (i = 0; i < nseq; i++)
	{
		int lenp = len[i] - oln1;
		pos[i] = rand() % lenp;
		ori[i] = rand() % 2;
	}
	for (i = 0; i < 16; i++)deg[i] = 0;
	i = 0;
	do
	{
		int r = rand() % 16;
		if (deg[r] == oln2)continue;
		deg[r]++;
		i++;
	} while (i < size);
	int t = 0;
	for (i = 0; i < 16; i++)
	{
		if (deg[i] > 0)
		{
			int take_pos[POPSIZE];
			for (j = 0; j < deg[i]; j++)take_pos[j] = 1;
			for (j = deg[i]; j < oln2; j++)take_pos[j] = 0;
			//		printf("DO\t");
			//	for(j=0;j<oln1;j++)printf("%d\t",take_pos[j]);printf("\n");
			BigMixI(take_pos, oln2);
			//	printf("PO\t");
			//	for(j=0;j<oln1;j++)printf("%d\t",take_pos[j]);printf("\n");
			for (j = 0; j < oln2; j++)
			{
				if (take_pos[j] == 1)
				{
					tot[t].num = i;
					tot[t].sta = tot[t].end = j;
					t++;
				}
			}
			t -= deg[i];
			int jend = deg[i] - 1;
			for (j = 0; j < deg[i]; j++)
			{
				int rr = rand() % 2;
				int rlen = reg_max - tot[t].end + tot[t].sta;
				if (rr == 0)
				{
					int prev_pos;
					if (j == 0)prev_pos = -1;
					else prev_pos = tot[t - 1].end;
					int spac = tot[t].sta - prev_pos;
					if (spac > 1)
					{
						if (spac > rlen)spac = rlen;
						int sh = rand() % spac;
						tot[t].sta -= sh;
					}
				}
				else
				{
					int next_pos;
					if (j != jend)next_pos = tot[t + 1].sta;
					else next_pos = oln1;
					int spac = next_pos - tot[t].end;
					if (spac > 1)
					{
						if (spac > rlen)spac = rlen;
						int sh = rand() % spac;
						tot[t].end += sh;
					}
				}
				odg[tot[t].end - tot[t].sta]++;
				t++;
			}
		}
	}
}

int town::mem_in(int nseq)
{
	pos = new int[nseq];
	if (pos == NULL) return -1;
	ori = new int[nseq];
	if (ori == NULL) return -1;
	return 1;
}
void town::mem_out(void)
{
	delete[] pos;
	delete[] ori;
}
int town::sum(int j)
{
	int i;
	int ret = 0;
	if (j == 0) return 0;
	for (i = 0; i<j; i++)ret += deg[i];
	return ret;
}
void town::fprint_all(char *file, char *add)
{
	int i;
	FILE *out;
	char file_out[500];
	strcpy(file_out, file);
	strcat(file_out, add);
	if ((out = fopen(file_out, "at")) == NULL)
	{
		printf("Ouput file can't be opened!\n");
		exit(1);
	}
	fprintf(out, "\t//%s\n", file);
	fprintf(out, "%d\t%f\t", size, fit);
	fprintf(out, "\n");
	for (i = 0; i<size; i++)
	{
		if (i == 0 || (i != 0 && (tot[i].num != tot[i - 1].num)))
		{
			fprintf(out, "%d\t", deg[tot[i].num]);
		}
		fprintf(out, "[%d;%d]%s\t", tot[i].sta, tot[i].end, s[tot[i].num].oli);
		if (tot[i].num != tot[i + 1].num)fprintf(out, "\n");
	}
	fclose(out);
}
void town::fprint_allfi(char *file, char *add, int len, double c0, double *buf)
{
	int i;
	FILE *out;
	char file_out[500];
	strcpy(file_out, file);
	strcat(file_out, add);
	if ((out = fopen(file_out, "at")) == NULL)
	{
		printf("Ouput file can't be opened!\n");
		exit(1);
	}
	fprintf(out, "if(strcmp(site,\"%s\")==0)\n", file);// 0 - one window system
	fprintf(out, "{\n\tcity a1 = {");
	fprintf(out, "\n\t//site,size,len,c0\n\t\"%s\",%d,%d,%f,0.05,\n\t", file, size, len, c0);
	fprintf(out, "{//buf,sta,end,num\t");
	for (i = 0; odg[i] != -1; i++)fprintf(out, "%d ", odg[i]); fprintf(out, "\n\t");
	for (i = 0; i<size; i++)
	{
		fprintf(out, "{%9f,%d,%d,%d}", buf[i], tot[i].sta, tot[i].end, tot[i].num);
		if (i == size - 1)
		{
			fprintf(out, "}};\t//%s\n", s[tot[i].num].oli);
			break;
		}
		else fprintf(out, ",");
		if (tot[i].num != tot[i + 1].num)fprintf(out, "\t//%s\n\t", s[tot[i].num].oli);
	}
	fprintf(out, "\ta1.get_copy(&sta);\n}\n");
	fclose(out);
}
void town::fprint_allfi_mat(char *file, char *add, char *name, int len, double c0, double *buf)
{
	int i;
	FILE *out;
	char file_out[500];
	strcpy(file_out, file);
	strcat(file_out, add);

	if ((out = fopen(file_out, "at")) == NULL)
	{
		printf("Ouput file can't be opened!\n");
		exit(1);
	}
	fprintf(out, "%s\n", name);// 0 - one window system	
	fprintf(out, "%d\tLPD count\n", size);
	fprintf(out, "%d\tModel length\n", len);
	fprintf(out, "%f\tCoefficient\n", c0);
	for (i = 0; i<size; i++)
	{
		fprintf(out, "%d\t%d\t%9f\t", tot[i].sta, tot[i].end, buf[i]);
		fprintf(out, "%d\t%s\n", tot[i].num, s[tot[i].num].oli);
	}
	fclose(out);
}
void town::fprint_seq(char *file, int len, int nseq, char ***seq, double *best_sco)
{
	int i, j, k;
	char d[POPSIZE], dir[] = "+-";
	FILE *out;
	if ((out = fopen(file, "at")) == NULL)
	{
		printf("Ouput file can't be opened!\n");
		exit(1);
	}
	for (i = 0; i<nseq; i++)
	{
		k = 0;
		int x1 = pos[i], x2 = pos[i] + len - 1;
		if (ori[i] == 0)
		{
			for (j = x1; j <= x2; j++)d[k++] = seq[0][i][j];
		}
		else
		{
			for (j = x2; j >= x1; j--)d[k++] = seq[1][i][j];
		}
		d[len] = '\0';
		fprintf(out, "%d\t%d\t%d\t%c\t%f\t%s\n", i + 1, x1, x2, dir[ori[i]], best_sco[i], d);
	}
	fclose(out);
}
void GetWords(int word, int size0, int size, char *w0)
{
	int i, j, zna, che, otn;
	for (i = size0; i<size0 + size; i++)
	{
		che = i - size0;
		zna = size / 4;
		memset(s[i].oli, 0, sizeof(s[i].oli));
		for (j = 0; j<word; j++)
		{
			otn = che / zna;
			s[i].oli[j] = w0[otn];
			che -= otn*zna;
			zna /= 4;
		}
		s[i].oli[word] = '\0';
		//printf("%d\t%s\n",i,s[i].oli);
	}
}
/*
void GetSost2(char *d, stru5 a)
{
int l, i, j, k;
int t[5]={10,11,20,21,22};
int nt=2;
for(i=0;i<a.size;i++)sost[i]=0;
int len=strlen(d), word;
double x, dx;
int left0, right, left;
*/
void DelHole(char *str)
{
	char *hole;
	hole = strstr(str, "\n");
	if (hole != NULL) *hole = 0;
}
char *TransStr(char *d)
{
	int i, c, lens;
	lens = strlen(d);
	for (i = 0; i<lens; i++)
	{
		c = int(d[i]);
		if (c<97) d[i] = char(c + 32);
		//else break;
	}
	return(d);
}
int CheckStr(char *file, char *d, int n, int print)
{
	int i, len, ret;
	len = strlen(d);
	ret = 1;
	for (i = 0; i < len; i++)
	{
		if (strchr("atgcATGC\n", (int)d[i]) != NULL)continue;
		if (print == 1)printf("File %s; sequence %d position %d (%c) bad. Sequence deleted!\n", file, n, i + 1, d[i]);
		ret = -1;
		break;
	}
	return(ret);
}
int IdeLet(char c)
{
	int ret;
	switch (c){
	case 'a': ret = 0; break;
	case 'c': ret = 1; break;
	case 'g': ret = 2; break;
	case 't': ret = 3; break;
	default: ret = -1;
	}
	return(ret);
}
void MinusStr1(int size, int j, int j0, double buf, double b[POPSIZE][POPSIZE])
{
	int k;
	b[j][j0] = 0;
	for (k = j0 + 1; k<size; k++)
	{
		if (b[j0][k] == 0)continue;
		b[j][k] -= buf*b[j0][k];
	}
}
void MinusStr2(int size, int j, int j0, double buf, double b[POPSIZE][POPSIZE])
{
	int k;
	b[j][j0] = 0;
	for (k = j - 1; k >= j0; k--)
	{
		if (b[j0][k] == 0)continue;
		b[j][k] -= buf*b[j0][k];
	}
}
void MinusStr(int size, int j, int j0, double buf, double b[POPSIZE][POPSIZE])
{
	int k;
	for (k = 0; k<size; k++)
	{
		if (b[j0][k] == 0)continue;
		b[j][k] -= buf*b[j0][k];
	}
}
int BackMat(int size)
{
	int i, j;
	double buf, b[POPSIZE][POPSIZE];
	if (size == 1)
	{
		uw[0][0] = 1 / uw[0][0];
		return 1;
	}
	for (j = 0; j<size; j++)
	{
		for (i = 0; i<size; i++)
		{
			b[i][j] = 0;
		}
	}
	for (i = 0; i<size; i++)b[i][i] = 1;
	for (j = 0; j<size - 1; j++)
	{
		for (i = j + 1; i<size; i++)
		{
			if (uw[i][j] == 0)continue;
			buf = uw[i][j] / uw[j][j];
			MinusStr1(size, i, j, buf, uw);
			MinusStr(size, i, j, buf, b);
			if (fabs(buf)<1e-035)
			{
				printf("\n1back %g(%d,%d) bij%g(%d,%d)\t", uw[i][j], i, j, b[i][j], i, j);
				//printf("%g(%d)\n",uw[j][j],j);
				return -1;
				//
			}
			uw[i][j] = 0;
		}
	}
	for (j = size - 1; j>0; j--)
	{
		for (i = j - 1; i >= 0; i--)
		{
			if (uw[i][j] == 0)continue;
			buf = uw[i][j] / uw[j][j];
			MinusStr2(size, i, j, buf, uw);
			MinusStr(size, i, j, buf, b);
			if (fabs(buf)<1e-035)
			{
				printf("\n2back %g(%d,%d) bij%g(%d,%d)\t", uw[i][j], i, j, b[i][j], i, j);
				//printf("%g(%d)\n",uw[j][j],j);
				//
				return -1;
			}
			uw[i][j] = 0;
		}
	}
	for (i = 0; i<size; i++)
	{
		buf = uw[i][i];
		if (fabs(buf)<1e-030)
		{
			printf("\nback3 %g %d\t", uw[i][i], i);
			return -1;
		}
		for (j = 0; j<size; j++)
		{
			uw[i][j] /= buf;
			b[i][j] /= buf;
		}
	}
	/*
	for(i=0;i<size;i++)
	{

	for(j=0;j<size;j++)
	{
	printf("%.2e ",uw[i][j]);
	}
	printf("\n");
	}
	printf("\n");
	for(i=0;i<size;i++)
	{
	for(j=0;j<size;j++)
	{
	printf("%.1e ",b[i][j]);
	}
	printf("\n");
	}
	printf("\n");
	int n;
	for(i=0;i<size;i++)
	{
	for(j=0;j<size;j++)
	{
	double x=0;
	for(n=0;n<size;n++)x+=uw0[j][n]*b[n][i];
	printf("%.1e ",x);
	}
	printf("\n");
	}
	*/
	for (i = 0; i<size; i++)
	{
		for (j = 0; j<size; j++)
		{
			uw[i][j] = b[i][j];
		}
	}
	return 1;
}
/*
void Test(char **peak_real, int *len, int num, int tryx)
{
char d[SEQLEN];
int k;
printf("Test %d\tlen=%d\n",tryx,len[num]);
for(k=0;k<len[num];k++)
{
d[k]=peak_real[num][k];
printf("%c",d[k]);
}
printf("\n");
}
*/
double EvalMahFIT(town *a, int n_train, int ***seq, double **dav, double **dcv, double **frp)
{
	int k, n, m, b;
	double f1[POPSIZE];//av[POPSIZE] buf[POPSIZE],
	double df[POPSIZE];

	for (k = 0; k < a->size; k++)
	{
		for (n = 0; n < a->size; n++)uw[k][n] = 0;
	}
	for (k = 0; k < a->size; k++)df[k] = 0;
	for (m = 0; m< n_train; m++)
	{		
		int ori = a->ori[m];
		int pos = a->pos[m];
		double fs[POPSIZE];
		for (k = 0; k < a->size; k++)
		{
			int rlenk = (a->tot[k].end - a->tot[k].sta + 1);
			fs[k] = 0;
			for (n = a->tot[k].sta; n <= a->tot[k].end; n++)
			{
				if (a->tot[k].num == seq[ori][m][n + pos])fs[k]++;
			}
			fs[k] /= rlenk;
			frp[m][k] = fs[k];
		}
		for (k = 0; k < a->size; k++)
		{
			for (n = 0; n < a->size; n++)
			{
				uw[k][n] += fs[k] * fs[n];
			}
			df[k] += fs[k];
		}
	}
	for (k = 0; k < a->size; k++)for (n = 0; n < a->size; n++)uw[k][n] /= n_train;
	for (k = 0; k < a->size; k++)df[k] /= n_train;
	for (k = 0; k < a->size; k++)f1[k] = df[k];
	for (k = 0; k < a->size; k++)
	{
		for (n = 0; n < a->size; n++)
		{
			uw[k][n] -= df[k] * df[n];
		}
	}
	for (k = 0; k < a->size; k++)
	{
		uw[k][k] += dcv[a->tot[k].end - a->tot[k].sta][a->tot[k].num];
	}
	for (k = 0; k < a->size; k++)for (n = 0; n < a->size; n++)uw[k][n] /= 2;
	for (k = 0; k < a->size; k++)
	{
		df[k] = df[k] - dav[a->tot[k].end - a->tot[k].sta][a->tot[k].num];
	}
	if (BackMat(a->size) == -1) { a->fit = 0; return 0; }
	a->mah = 0;
	for (k = 0; k < a->size; k++)
	{
		double buf = 0;
		for (n = 0; n < a->size; n++)buf += uw[k][n] * df[n];
		a->mah += buf * df[k];
	}
	//a->ave = a->std = 1;	
	double xtr = 0, xtr2 = 0;
	for (b = 0; b < n_train; b++)
	{
		/*m = xporti[b];
		int ori = a->ori[m];
		int pos = a->pos[m];
		double fs[POPSIZE];
		for (k = 0; k < a->size; k++)
		{
		int rlenk = (a->tot[k].end - a->tot[k].sta + 1);
		fs[k] = 0;
		for (n = a->tot[k].sta; n <= a->tot[k].end; n++)
		{
		if (a->tot[k].num == seq[ori][m][n + pos])fs[k]++;
		}
		fs[k] /= rlenk;
		}*/
		double mah1 = 0;
		for (k = 0; k < a->size; k++)
		{
			double buf1 = 0;
			for (n = 0; n < a->size; n++)buf1 += uw[k][n] * (frp[b][n] - f1[n]);
			mah1 += buf1 * (frp[b][k] - f1[k]);
		}
		xtr += mah1;
		xtr2 += mah1*mah1;
	}
	xtr /= n_train;
	xtr2 /= n_train;
	a->ave = xtr;
	a->std = xtr2 - xtr*xtr;
	if (a->std <= 0.001 || a->std > 1000)
	{
		a->fit = 0;
		return 0;
	}
	a->std = sqrt(a->std);
	a->fit = a->mah * a->ave / a->std;
	frp = NULL;
	return a->fit;
}
double EvalMahFITTrain(town *a, int nseq, int ***seq, char *file, int olen, int *len, char ***peak_real, town_ext *best_sel_ext, double *best_sco, double **dav, double **dcv)
{
	int k, n, m;
	double ret = 0;
	double av[POPSIZE], buf[POPSIZE];
	for (k = 0; k<a->size; k++)
	{
		uw0[k] = dcv[a->tot[k].end - a->tot[k].sta][a->tot[k].num];
		//printf("%f\t[%d;%d]%s\t[%d;%d]%s\t%d\t%d\t%d\n",uw0[k][n],a->tot[k].sta,a->tot[k].end,s[a->tot[k].num].oli,a->tot[n].sta,a->tot[n].end,s[a->tot[n].num].oli,k1+2,k2+1,k3);			
	}
	for (k = 0; k<a->size; k++)
	{
		for (n = 0; n<a->size; n++)uw[k][n] = 0;
	}
	for (k = 0; k<a->size; k++)df[k] = 0;
	for (m = 0; m<nseq; m++)
	{
		int ori = a->ori[m];
		int pos = a->pos[m];
		double fs[POPSIZE];
		for (k = 0; k<a->size; k++)
		{
			int rlenk = (a->tot[k].end - a->tot[k].sta + 1);
			fs[k] = 0;
			for (n = a->tot[k].sta; n <= a->tot[k].end; n++)
			{
				if (a->tot[k].num == seq[ori][m][n + pos])fs[k]++;
			}
			fs[k] /= rlenk;
		}
		for (k = 0; k<a->size; k++)
		{
			for (n = 0; n<a->size; n++)
			{
				uw[k][n] += fs[k] * fs[n];
			}
			df[k] += fs[k];
		}
	}
	for (k = 0; k<a->size; k++)for (n = 0; n<a->size; n++)uw[k][n] /= nseq;
	for (k = 0; k<a->size; k++)df[k] /= nseq;
	for (k = 0; k<a->size; k++)
	{
		for (n = 0; n<a->size; n++)
		{
			uw[k][n] -= df[k] * df[n];
		}
	}
	for (k = 0; k<a->size; k++)
	{
		uw[k][k] += uw0[k];
	}
	for (k = 0; k<a->size; k++)for (n = 0; n<a->size; n++)uw[k][n] /= 2;
	{
		for (k = 0; k<a->size; k++)
		{
			av[k] = (df[k] + dav[a->tot[k].end - a->tot[k].sta][a->tot[k].num]) / 2;
			df[k] = df[k] - dav[a->tot[k].end - a->tot[k].sta][a->tot[k].num];
		}
	}
	if (BackMat(a->size) == -1){ a->fit = 0; return 0; }
	a->fit = 0;
	//double euc_buf=0;
	for (k = 0; k<a->size; k++)
	{
		buf[k] = 0;
		for (n = 0; n<a->size; n++)buf[k] += uw[k][n] * df[n];
		a->fit += buf[k] * df[k];
		//	euc_buf+=buf[k]*buf[k];
		//printf("%f\t%f\t%f\n",a->fit,buf,df[k]);
	}
	//	euc_buf=sqrt(euc_buf);
	//a->fit/=euc_buf;
	a->fit /= 2;
	double c0 = 0;
	for (k = 0; k<a->size; k++)
	{
		buf[k] /= a->fit;
		c0 -= av[k] * buf[k];
	}
	char ext2[] = "_end2";
	a->fprint_allfi(file, ext2, olen, c0, buf);
	best_sel_ext->get_copy(c0, buf, a->size);
	int b, o;
	for (b = 0; b<nseq; b++)
	{
		double fs[POPSIZE];
		int lenp = len[b] - olen + 1;
		double sco_pos = -1000;
		for (m = 0; m<lenp; m++)
		{
			for (o = 0; o<2; o++)
			{
				double sco = c0;
				for (k = 0; k<a->size; k++)
				{
					fs[k] = 0;
					int rlenk = (a->tot[k].end - a->tot[k].sta + 1);
					for (n = a->tot[k].sta; n <= a->tot[k].end; n++)
					{
						if (a->tot[k].num == seq[o][b][n + m])fs[k]++;
					}
					if (fs[k] != 0)
					{
						fs[k] /= rlenk;
						sco += buf[k] * fs[k];
					}
				}
				sco = 1 - fabs(1 - sco);
				if (m == 0 && o == 0)sco_pos = sco;
				else
				{
					if (sco > sco_pos)
					{
						sco_pos = sco;
					}
				}
			}
		}
		best_sco[b] = sco_pos;
	}
	return a->fit;
}
int MutOlig0(town *a)
{
	int cy = 0, r1, num1, i, gom;
	do
	{
		gom = 0;
		r1 = rand() % a->size;
		num1 = rand() % 15;
		if (num1 >= a->tot[r1].num)num1++;
		int i0 = a->sum(num1);
		for (i = i0; i<i0 + a->deg[num1]; i++)
		{
			//if(a->tot[i].num!=num1)continue;
			//if(i==r1)continue;
			if (a->tot[r1].sta >= a->tot[i].sta && a->tot[r1].sta <= a->tot[i].end){ gom = 1; break; }
			if (a->tot[r1].end >= a->tot[i].sta && a->tot[r1].end <= a->tot[i].end){ gom = 1; break; }
			if (a->tot[i].sta >= a->tot[r1].sta && a->tot[i].sta <= a->tot[r1].end){ gom = 1; break; }
			if (a->tot[i].end >= a->tot[r1].sta && a->tot[i].end <= a->tot[r1].end){ gom = 1; break; }
		}
		cy++;
		if (cy>10)return -1;
	} while (gom == 1);
	if (gom == 1)return -1;
	a->deg[a->tot[r1].num]--;
	a->deg[num1]++;
	a->tot[r1].num = num1;
	return r1;
}
int CTown(town *a, int *n, int len, int olen)
{
	int j, k;
	n[0] = n[1] = -1;
	int i, gom, ns = 0;
	int odgs[] = { 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0,
		0, 0, 0, 0, 0,
		0, 0, 0, 0, 0,
		0, 0, 0, 0, 0,
		0, 0, 0, 0, 0 };
	len -= 1;
	for (j = 0; j<a->size; j++)
	{
		if (a->tot[j].sta>a->tot[j].end)
		{
			printf("LEN STA %d\tEND %d\n", a->tot[j].sta, a->tot[j].end);
			return -1;
		}
		if (a->tot[j].end - a->tot[j].sta + 1>olen)
		{
			printf("LEN STA %d\tEND %d\n", a->tot[j].sta, a->tot[j].end);
			return -1;
		}
	}
	for (j = 0; j<a->size; j++)
	{

		if (a->tot[j].sta >= len || a->tot[j].sta<0){ n[0] = j; return -1; }
		if (a->tot[j].end >= len || a->tot[j].end<0){ n[0] = j; return -1; }
		odgs[a->tot[j].end - a->tot[j].sta]++;
	}
	for (i = 0; i<olen; i++)
	{
		if (a->odg[i] != odgs[i])
		{
			printf("Oli %d Announced %d Really %d", i + 1, a->odg[i], odgs[i]);
			return -1;
		}
	}
	for (j = 0; j<16; j++)
	{
		if (a->deg[j] == 0)continue;
		//if(a->deg[j]==1){ns++;continue;}
		gom = 0;
		for (i = ns; i<ns + a->deg[j]; i++)
		{
			for (k = ns; k<ns + a->deg[j]; k++)
			{
				if (i == k)continue;
				if (a->tot[i].sta >= a->tot[k].sta && a->tot[i].sta <= a->tot[k].end)
				{
					gom = 1; n[0] = i; n[1] = k; break;
				}
				if (a->tot[i].end >= a->tot[k].sta && a->tot[i].end <= a->tot[k].end)
				{
					gom = 1; n[0] = i; n[1] = k; break;
				}
			}
			if (gom == 1)break;
		}
		if (gom == 1)
		{
			printf("Overlap %d\t%d\n", n[0], n[1]);
			a->tot[n[0]].print_all();
			a->tot[n[1]].print_all();
			return -1;
		}
		else ns += a->deg[j];
	}
	int sum = 0;
	for (i = 0; i<olen; i++)
	{
		if (a->odg[i]<0)
		{
			printf("Oli %d\todg = %d\n", i, a->odg[i]);
			return -1;
		}
		sum += a->odg[i];
	}
	if (sum != a->size)
	{
		printf("Sum by odg %d\tsize = %d\n", sum, a->size);
		return -1;
	}
	return 1;
}
int MutCry0(town *a, int oln, int reg_max)
{
	int n1;
	int a1 = 0, b1 = oln - 1, oln2 = oln - 2;
	int fl = 1, cy = 0;
	do
	{
		n1 = rand() % a->size;
		int n1m1 = n1 - 1;
		if (n1 != 0 && a->tot[n1m1].num == a->tot[n1].num)
		{
			if (a->tot[n1m1].end<oln2)a1 = a->tot[n1m1].end + 1;
			else
			{
				fl = -1;
				cy++;
				continue;
			}
		}
		else a1 = 0;
		int n1p1 = n1 + 1;
		if (n1 != a->size - 1 && a->tot[n1p1].num == a->tot[n1].num)
		{
			if (a->tot[n1p1].sta >= 1)b1 = a->tot[n1p1].sta - 1;
			else
			{
				fl = -1;
				cy++;
				continue;
			}
		}
		else
		{
			b1 = oln2;
		}
		if (cy>10)return -1;
	} while (fl == -1);
	if (a1 >= b1) return -1;
	int a2, b2, ab = Min(reg_max, b1 - a1 + 1);//max available len oliga 
	int w = rand() % ab;
	int loc = b1 - a1 - w + 1;
	if (loc <= 0)return -1;
	int da2 = rand() % loc;
	a2 = a1 + da2;
	b2 = a2 + w;
	if (a->tot[n1].sta == a2 && a->tot[n1].end == b2)return -1;
	int w0 = a->tot[n1].end - a->tot[n1].sta;
	//	printf("%d%d ",a->odg[w0],a->odg[w]);
	a->odg[w0]--;
	a->odg[w]++;
	//	printf("%d%d \n",a->odg[w0],a->odg[w]);
	a->tot[n1].sta = a2;
	a->tot[n1].end = b2;
	return n1;
}
int MutRegShift(town *a, int nseq, int *len, int olen, int &npeak, int &nori, int &npos)
{
	//printf("In Peak %d Ori %d Pos %d ",npeak,nori,npos);
	int r2, cy = 0;
	npeak = rand() % nseq;
	int lenr = len[npeak] - olen + 1;
	r2 = rand() % lenr;
	npos = r2;
	if (r2 == a->pos[npeak])
	{
		nori = a->ori[npeak] = 1 - a->ori[npeak];
		//printf("Out1 Peak %d Pos %d Ori %d",npeak, npos,nori);
		return 1;
	}
	else
	{
		a->pos[npeak] = r2;
		int r3 = rand() % 2;
		if (r3 == 1)nori = a->ori[npeak] = 1 - a->ori[npeak];
		else nori = a->ori[npeak];
		//printf("Out2 Peak %d Pos %d Ori %d",npeak, npos,nori);
		return 1;
	}
}int RecFeat(town a1, town a2, int(*cop)[2], int max)
{
	int ret = 0;
	int reg[2] = { a1.size / 16, a2.size / 16 };
	int cy, r1, r2, num1, num2, i, gom1, gom2, r0;
	do
	{
		cy = 0;
		do
		{
			gom1 = 0;
			r1 = rand() % (a1.size - ret);
			for (i = 0; i<ret; i++)
			{
				if (r1 >= cop[i][0])r1++;
			}
			num1 = a1.tot[r1].num;
			if (a2.deg[num1] == 0)break;
			int i0 = a2.sum(num1);
			for (i = i0; i<i0 + a2.deg[num1]; i++)
			{
				if (gom1>1)break;
				if (a2.tot[i].sta>a1.tot[r1].end)break;
				if (a2.tot[i].end<a1.tot[r1].sta)continue;
				if (a2.tot[i].sta >= a1.tot[r1].sta && a2.tot[i].sta <= a1.tot[r1].end){ gom1++; r0 = i; continue; }
				if (a2.tot[i].end >= a1.tot[r1].sta && a2.tot[i].end <= a1.tot[r1].end){ gom1++; r0 = i; continue; }
				if (a1.tot[r1].sta >= a2.tot[i].sta && a1.tot[r1].sta <= a2.tot[i].end){ gom1++; r0 = i; continue; }
				if (a1.tot[r1].end >= a2.tot[i].sta && a1.tot[r1].end <= a2.tot[i].end){ gom1++; r0 = i; continue; }
			}
			if (gom1 == 1)
			{
				if (a1.tot[r1].end == a2.tot[r0].end && a1.tot[r1].sta == a2.tot[r0].sta)
				{
					gom1++;
				}
			}
			if (cy>reg[0])return ret;
			cy++;
		} while (gom1>1);
		cy = 0;
		do
		{
			gom2 = 0;
			if (gom1 == 0)r2 = rand() % (a2.size - ret);
			else r2 = r0;
			for (i = 0; i<ret; i++)
			{
				if (r2 >= cop[i][1])r2++;
			}
			num2 = a2.tot[r2].num;
			if (a1.deg[num2] == 0)break;
			if (num1 == num2 && (a1.tot[r1].end == a2.tot[r2].end && a1.tot[r1].sta == a2.tot[r2].sta))continue;
			int i0 = a1.sum(num2);
			for (i = i0; i<i0 + a1.deg[num2]; i++)
			{
				if (gom2>1)break;
				if (a1.tot[i].sta>a2.tot[r2].end)break;
				if (a1.tot[i].end<a2.tot[r2].sta)continue;
				if (a1.tot[i].sta >= a2.tot[r2].sta && a1.tot[i].sta <= a2.tot[r2].end){ gom2++; continue; }
				if (a1.tot[i].end >= a2.tot[r2].sta && a1.tot[i].end <= a2.tot[r2].end){ gom2++; continue; }
				if (a2.tot[r2].sta >= a1.tot[i].sta && a2.tot[r2].sta <= a1.tot[i].end){ gom2++; continue; }
				if (a2.tot[r2].end >= a1.tot[i].sta && a2.tot[r2].end <= a1.tot[i].end){ gom2++; continue; }
			}
			if (gom1 == 1 && gom2 == 1)break;
			if (cy>reg[1])return ret;
			cy++;
		} while (gom2 != 0);
		cop[ret][0] = r1;
		cop[ret][1] = r2;
		//		a1.tot[r1].print_all();
		//		a2.tot[r2].print_all();
		//		printf("\n");
		ret++;
	} while (ret<max);
	//	a1.print_all();
	//	a2.print_all();
	return ret;
}
//original: simple exchange
int Reco2_Original(town *a1, town *a2, int *nsi, int *num)
{
	int cy = 0, r1, r2, num1, num2, i, gom1, gom2, r0;
	//int size1=Min(a1->size,a2->size);
	do
	{
		gom1 = 0;
		r1 = rand() % a1->size;
		num1 = a1->tot[r1].num;
		if (a2->deg[num1] == 0)break;
		int i0 = a2->sum(num1);
		for (i = i0; i<i0 + a2->deg[num1]; i++)
		{
			if (gom1>1)break;
			if (a2->tot[i].sta>a1->tot[r1].end)break;
			if (a2->tot[i].end<a1->tot[r1].sta)continue;
			if (a2->tot[i].sta >= a1->tot[r1].sta && a2->tot[i].sta <= a1->tot[r1].end){ gom1++; r0 = i; continue; }
			if (a2->tot[i].end >= a1->tot[r1].sta && a2->tot[i].end <= a1->tot[r1].end){ gom1++; r0 = i; continue; }
			if (a1->tot[r1].sta >= a2->tot[i].sta && a1->tot[r1].sta <= a2->tot[i].end){ gom1++; r0 = i; continue; }
			if (a1->tot[r1].end >= a2->tot[i].sta && a1->tot[r1].end <= a2->tot[i].end){ gom1++; r0 = i; continue; }
		}
		if (gom1 == 1)
		{
			if (a1->tot[r1].end == a2->tot[r0].end && a1->tot[r1].sta == a2->tot[r0].sta)
			{
				gom1++;
			}
		}
		if (cy>10)return -1;
		cy++;
	} while (gom1>1);
	cy = 0;
	do
	{
		gom2 = 0;
		if (gom1 == 0)r2 = rand() % a2->size;
		else r2 = r0;
		num2 = a2->tot[r2].num;
		if (a1->deg[num2] == 0)break;
		if (num1 == num2 && (a1->tot[r1].end == a2->tot[r2].end && a1->tot[r1].sta == a2->tot[r2].sta))continue;
		int i0 = a1->sum(num2);
		for (i = i0; i<i0 + a1->deg[num2]; i++)
		{
			if (gom2>1)break;
			if (a1->tot[i].sta>a2->tot[r2].end)break;
			if (a1->tot[i].end<a2->tot[r2].sta)continue;
			if (a1->tot[i].sta >= a2->tot[r2].sta && a1->tot[i].sta <= a2->tot[r2].end){ gom2++; continue; }
			if (a1->tot[i].end >= a2->tot[r2].sta && a1->tot[i].end <= a2->tot[r2].end){ gom2++; continue; }
			if (a2->tot[r2].sta >= a1->tot[i].sta && a2->tot[r2].sta <= a1->tot[i].end){ gom2++; continue; }
			if (a2->tot[r2].end >= a1->tot[i].sta && a2->tot[r2].end <= a1->tot[i].end){ gom2++; continue; }
		}
		if (gom1 == 1 && gom2 == 1)break;
		if (cy>10)return -1;
		cy++;
	} while (gom2 != 0);
	a1->deg[num1]--;
	a2->deg[num2]--;
	a2->deg[num1]++;
	a1->deg[num2]++;
	int w1 = a1->tot[r1].end - a1->tot[r1].sta;
	int w2 = a2->tot[r2].end - a2->tot[r2].sta;
	a1->odg[w1]--;
	a1->odg[w2]++;
	a2->odg[w2]--;
	a2->odg[w1]++;
	uno buf;
	a1->tot[r1].get_copy(&buf);
	a2->tot[r2].get_copy(&a1->tot[r1]);
	buf.get_copy(&a2->tot[r2]);
	num[0] = num1;
	num[1] = num2;
	nsi[0] = r1;
	nsi[1] = r2;
	return 1;
}
//economic: two slightly differnt feature exchanged
int Reco2_Economic(town *a1, town *a2, int *nsi, int *num)
{
	int i0, j, r1, r2, num1, num2, i, gom1 = 0, gom2 = 0, ret;
	int size1 = Min(a1->size, a2->size);
	int ord[POPSIZE];
	for (i = 0; i<size1; i++)ord[i] = i;
	BigMixI(ord, size1);
	size1 /= 2;
	for (j = 0; j<size1; j++)
	{
		r1 = ord[j];
		gom1 = 0;
		num1 = a1->tot[r1].num;
		if (a2->deg[num1] == 0)continue;
		i0 = a2->sum(num1);
		if (a2->deg[num1] == 1)r2 = i0;
		else
		{
			for (i = i0; i<i0 + a2->deg[num1]; i++)
			{
				if (a2->tot[i].end<a1->tot[r1].sta)continue;
				if (gom1>1)break;
				if (a2->tot[i].sta>a1->tot[r1].end)break;
				if (a2->tot[i].sta >= a1->tot[r1].sta && a2->tot[i].sta <= a1->tot[r1].end){ gom1++; r2 = i; continue; }//1s 2s 1e
				if (a2->tot[i].end >= a1->tot[r1].sta && a2->tot[i].end <= a1->tot[r1].end){ gom1++; r2 = i; continue; }//1s 2e 1e
				if (a1->tot[r1].sta >= a2->tot[i].sta && a1->tot[r1].sta <= a2->tot[i].end){ gom1++; r2 = i; continue; }//2s 1s 2e
				if (a1->tot[r1].end >= a2->tot[i].sta && a1->tot[r1].end <= a2->tot[i].end){ gom1++; r2 = i; continue; }//2s 1e 2e
			}
		}
		if (gom1>1)continue;
		if (gom1 == 0)
		{
			gom1++;
			r2 = i0 + rand() % a2->deg[num1];
		}
		if (a1->tot[r1].num != a2->tot[r2].num)continue;
		if (a1->tot[r1].end == a2->tot[r2].end && a1->tot[r1].sta == a2->tot[r2].sta)continue;
		/*		if(r2!=0 && a2->tot[r2-1].end>=a1->tot[r1].sta)continue;
		if(r1!=0 && a1->tot[r1-1].end>=a2->tot[r2].sta)continue;
		if(r1!=a1->size-1 && a1->tot[r1+1].sta<=a2->tot[r2].end)continue;
		if(r2!=a2->size-1 && a2->tot[r2+1].sta<=a1->tot[r1].end)continue;
		*/
		num2 = a2->tot[r2].num;
		gom2 = 1;
		i0 = a1->sum(num2);
		for (i = i0; i<i0 + a1->deg[num2]; i++)
		{
			if (i == r1)continue;
			if (gom2>1)break;
			if (a1->tot[i].sta>a2->tot[r2].end)break;
			if (a1->tot[i].end<a2->tot[r2].sta)continue;
			if (a1->tot[i].sta >= a2->tot[r2].sta && a1->tot[i].sta <= a2->tot[r2].end){ gom2++; continue; }
			if (a1->tot[i].end >= a2->tot[r2].sta && a1->tot[i].end <= a2->tot[r2].end){ gom2++; continue; }
			if (a2->tot[r2].sta >= a1->tot[i].sta && a2->tot[r2].sta <= a1->tot[i].end){ gom2++; continue; }
			if (a2->tot[r2].end >= a1->tot[i].sta && a2->tot[r2].end <= a1->tot[i].end){ gom2++; continue; }
		}
		ret = j;
		if (gom1 == 1 && gom2 == 1)break;
	}
	if (gom1 != 1 || gom2 != 1)return -1;
	if ((r1<0 || r1>a1->size - 1) || (r2<0 || r2>a2->size - 1))
	{
		puts("Error");
	}
	a1->deg[num1]--;
	a2->deg[num2]--;
	a2->deg[num1]++;
	a1->deg[num2]++;
	int w1 = a1->tot[r1].end - a1->tot[r1].sta;
	int w2 = a2->tot[r2].end - a2->tot[r2].sta;
	a1->odg[w1]--;
	a1->odg[w2]++;
	a2->odg[w2]--;
	a2->odg[w1]++;
	uno buf;
	a1->tot[r1].get_copy(&buf);
	a2->tot[r2].get_copy(&a1->tot[r1]);
	buf.get_copy(&a2->tot[r2]);
	num[0] = num1;
	num[1] = num2;
	nsi[0] = r1;
	nsi[1] = r2;
	return ret;
}
//total: exchange of subsets of one dinucl.type
int Reco2_One_dinucleotide_full(town *a1, town *a2)
{
	int gom = -1, i, j, w1, w2, i1, i2;
	int size1 = 16;
	int ord[16];
	for (i = 0; i<size1; i++)ord[i] = i;
	BigMixI(ord, size1);
	for (j = 0; j<size1; j++)
	{
		if (a1->deg[ord[j]] == 0)continue;
		if (a1->deg[ord[j]] != a2->deg[ord[j]])continue;
		i1 = a1->sum(ord[j]);
		i2 = a2->sum(ord[j]);
		for (i = 0; i<a1->deg[ord[j]]; i++)
		{
			int k1 = i1 + i;
			int k2 = i2 + i;
			if (a1->tot[k1].sta != a2->tot[k2].sta){ gom = ord[j]; break; }
			if (a1->tot[k1].end != a2->tot[k2].end){ gom = ord[j]; break; }
		}
		if (gom != -1)break;
	}
	if (gom == -1)return -1;
	uno buf;
	for (j = 0; j<a1->deg[gom]; j++)
	{
		int k1 = i1 + j;
		int k2 = i2 + j;
		w1 = a1->tot[k1].end - a1->tot[k1].sta;
		w2 = a2->tot[k2].end - a2->tot[k2].sta;
		a1->odg[w1]--;
		a1->odg[w2]++;
		a2->odg[w2]--;
		a2->odg[w1]++;
		a1->tot[k1].get_copy(&buf);
		a2->tot[k2].get_copy(&a1->tot[k1]);
		buf.get_copy(&a2->tot[k2]);
	}
	return 1;
}
//canonical: : classic recombination within one dinucl.type
int Reco2_One_dinucleotide_local(town *a1, town *a2, int olen, int reg_max)
{
	int i, j, j0, w1, w2, i1, i2;
	int size1 = 16;
	int ord[16];
	int rpt = 0, rp[SEQLEN];
	for (i = 0; i<size1; i++)ord[i] = i;
	BigMixI(ord, size1);
	for (j = 0; j<size1; j++)
	{
		if (a1->deg[ord[j]] <= 1)continue;
		if (a1->deg[ord[j]] != a2->deg[ord[j]])continue;
		i1 = a1->sum(ord[j]);
		i2 = a2->sum(ord[j]);
		for (i = 0; i<a1->deg[ord[j]]; i++)
		{
			int k1 = i1 + i;
			int k2 = i2 + i;
			if (a1->tot[k1].sta != a2->tot[k2].sta || a1->tot[k1].end != a2->tot[k2].end)
			{
				if ((a1->tot[k1].sta >= a2->tot[k2].sta && a1->tot[k1].sta <= a2->tot[k2].end) || (a1->tot[k1].end >= a2->tot[k2].sta && a1->tot[k1].end <= a2->tot[k2].end))
				{
					int sta1 = Max(a1->tot[k1].sta, a2->tot[k2].sta);
					int end1 = Min(a1->tot[k1].end, a2->tot[k2].end);
					int over = end1 - sta1;
					if (over<0)continue;
					int win = Max(a2->tot[k2].end - a1->tot[k1].sta, a1->tot[k1].end - a2->tot[k2].sta);
					if (win >= olen)continue;
					j0 = ord[j];
					rp[rpt++] = i;
				}
			}
		}
		if (rpt != 0)break;
	}
	if (rpt == 0)return -1;
	int rpt0 = rand() % rpt;
	{
		int k1 = i1 + rp[rpt0];
		int k2 = i2 + rp[rpt0];
		w1 = a1->tot[k1].end - a1->tot[k1].sta;
		w2 = a2->tot[k2].end - a2->tot[k2].sta;
		a1->odg[w1]--;
		a2->odg[w2]--;
		w1 = a2->tot[k2].end - a1->tot[k1].sta;
		if (w1<0 || w1 >= reg_max)return -1;
		w2 = a1->tot[k1].end - a2->tot[k2].sta;
		if (w2<0 || w2 >= reg_max)return -1;
		a1->odg[w1]++;
		a2->odg[w2]++;
		int end1 = a1->tot[k1].end;
		a1->tot[k1].end = a2->tot[k2].end;
		a2->tot[k2].end = end1;
	}
	uno buf;
	int cha = 0;
	for (j = rp[rpt0] + 1; j<a1->deg[j0]; j++)
	{
		int k1 = i1 + j;
		int k2 = i2 + j;
		w1 = a1->tot[k1].end - a1->tot[k1].sta;
		w2 = a2->tot[k2].end - a2->tot[k2].sta;
		a1->odg[w1]--;
		a1->odg[w2]++;
		a2->odg[w2]--;
		a2->odg[w1]++;
		a1->tot[k1].get_copy(&buf);
		a2->tot[k2].get_copy(&a1->tot[k1]);
		buf.get_copy(&a2->tot[k2]);
	}
	return 1;
}
int Reco2Peak(town *a1, town *a2, int n_train)
{
	int r1, r2 = 0;
	r1 = rand() % n_train;
	do
	{
		r2 = rand() % n_train;
	} while (r1 == r2);
	if (a1->pos[r1] == a2->pos[r2] && a1->ori[r1] == a2->ori[r2])return -1;
	MixI(&a1->ori[r1], &a2->ori[r1]);
	MixI(&a1->pos[r1], &a2->pos[r1]);
	MixI(&a1->ori[r2], &a2->ori[r2]);
	MixI(&a1->pos[r2], &a2->pos[r2]);
	return 1;
}

void MixPop(town *a, town *b)
{
	town c = *a;
	*a = *b;
	*b = c;
}
int GomTown(town a, town b, int nseq)
{
	int i;
	for (i = 0; i<a.size; i++)
	{
		if (b.tot[i].sta != a.tot[i].sta)return 0;
		if (b.tot[i].end != a.tot[i].end)return 0;
		if (b.tot[i].num != a.tot[i].num)return 0;
	}
	for (i = 0; i<nseq; i++)
	{
		if (b.pos[i] != a.pos[i])return 0;
		if (b.ori[i] != a.ori[i])return 0;
	}
	return -1;
}
void DelChar(char *str, char c)
{
	int i, lens, size;

	size = 0;
	lens = strlen(str);
	for (i = 0; i<lens; i++)
	{
		if (str[i] != c)str[size++] = str[i];
	}
	str[size] = '\0';
}
int Fun(char *d, town sta, int len0, double *p, double *rec_buf, double rec_c)
{
	int i, j, ret, len;
	int k;
	char d1[200];

	len = strlen(d);
	ret = 1;
	for (k = 0; k <= len - len0; k++)
	{
		for (j = 0; j<len0; j++)d1[j] = d[k + j]; d1[len0] = '\0';
		p[k] = 0;
		double sco = rec_c;
		for (j = 0; j<sta.size; j++)
		{
			int rlenj = (sta.tot[j].end - sta.tot[j].sta + 1);
			double fm = 0;
			for (i = sta.tot[j].sta; i <= sta.tot[j].end; i++)
			{
				int cod = 4 * IdeLet(d1[i]) + IdeLet(d1[i + 1]);
				if (sta.tot[j].num == cod){ fm++; break; }
			}
			if (fm != 0)
			{
				fm /= rlenj;
				sco += rec_buf[j] * fm;
			}
		}
		p[k] = 1 - fabs(sco - 1);
	}
	return 1;
}
int ComplStr(char *d)
{
	int i, len;
	len = strlen(d);
	char d1[SEQLEN];
	strcpy(d1, d);
	//	memset(d,0,sizeof(d));
	for (i = 0; i<len; i++)
	{
		switch (d1[len - i - 1])
		{
		case 'a':{d[i] = 't'; break; }
		case 't':{d[i] = 'a'; break; }
		case 'c':{d[i] = 'g'; break; }
		case 'g':{d[i] = 'c'; break; }
		default: d[i] = 'n';
		}
	}
	return 1;
}
void GetSost(char *d, int word, int *sost)
{
	int i, j, k, i_sost, let;
	char letter[5] = "acgt";
	int ten[6] = { 1, 4, 16, 64, 256, 1024 };
	int lens = strlen(d);
	int size = 1;
	for (k = 0; k<word; k++)size *= 4;
	for (i = 0; i<size; i++)sost[i] = 0;
	for (i = 0; i<lens - word + 1; i++)
	{
		i_sost = 0;
		for (j = word - 1; j >= 0; j--)
		{
			for (k = 0; k<4; k++)
			{
				if (d[i + j] == letter[k]){ let = k; break; }
			}
			i_sost += ten[word - 1 - j] * let;
		}
		sost[i_sost]++;
	}
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
	while (n >= 0)
	{
		if (fgets(l, sizeof(l), in) == NULL) fl = -1;
		if (*l == '\n' && fl != -1)continue;
		if (((*l == symbol) || (fl == -1)) && (fl != 0))
		{
			int lenx = strlen(d);
			int check = CheckStr(file, d, n, 1);
			if (lenx >= olen && check == 1)nseq++;
			if (fl == -1)
			{
				fclose(in);
				break;
			}
		}
		if (*l == symbol)
		{
			memset(head, 0, sizeof(head));
			DelHole(l);
			strcpy(head, l);
			fl = 0; continue;
		}
		if (fl == 0)
		{
			memset(d, 0, sizeof(d));
			DelHole(l);
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
		DelHole(l);
		strcat(d, l);
	}
}
void EvalLen(char *file, int *len, int olen)
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
	int nn = 0, n = 0;
	while (n >= 0)
	{
		if (fgets(l, sizeof(l), in) == NULL) fl = -1;
		if (*l == '\n' && fl != -1)continue;
		if (((*l == symbol) || (fl == -1)) && (fl != 0))
		{
			int lenx = strlen(d);
			int check = CheckStr(file, d, n, 0);
			if (lenx >= olen && check == 1)len[n++] = lenx;
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
			DelHole(l);
			strcpy(head, l);
			fl = 0; continue;
		}
		if (fl == 0)
		{
			memset(d, 0, sizeof(d));
			DelHole(l);
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
		DelHole(l);
		strcat(d, l);
	}
}
void ReadSeq(char *file, int nseq, int *len, int ***seq_real, char ***peak_real, int olen)
{
	char l[SEQLEN], d[2][SEQLEN], head[400];
	int fl = 0, i, j;
	FILE  *in;

	if ((in = fopen(file, "rt")) == NULL)
	{
		printf("ReadSeq! Input file %s can't be opened!\n", file);
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
			int lenx = strlen(d[0]);
			int check = CheckStr(file, d[0], nn, 0);
			nn++;
			if (lenx >= olen && check == 1)
			{
				TransStr(d[0]);
				d[0][len[n]] = '\0';
				strcpy(d[1], d[0]);
				ComplStr(d[1]);
				d[1][len[n]] = '\0';
				for (j = 0; j < 2; j++)
				{
					strcpy(peak_real[j][n], d[j]);
					peak_real[j][n][len[n]] = '\0';
					int cod[2];
					cod[0] = IdeLet(d[j][0]);
					int len1 = len[n] - 1;
					for (i = 0; i < len1; i++)
					{
						cod[1] = IdeLet(d[j][i + 1]);
						if (cod[0] != -1 && cod[1] != -1)seq_real[j][n][i] = 4 * cod[0] + cod[1];
						else seq_real[j][n][i] = -1;
						cod[0] = cod[1];
					}
				}
				n++;
			}
			else
			{
				if (lenx < olen)printf("peak %d (Len %d) ignored\n", n + 1, lenx);
				if (check == -1)printf("Unusual symbol, peak %d ignored\n%s\n", n + 1, d[0]);
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
			DelHole(l);
			strcpy(head, l);
			fl = 0; continue;
		}
		if (fl == 0)
		{
			memset(d[0], 0, sizeof(d[0]));
			DelHole(l);
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
		DelHole(l);
		strcat(d[0], l);
	}
}
void ReadSeqBack(char *file, int nseq, int *len, int olen, int reg_max, double **dav, double **dcv)
{
	char l[SEQLEN], d[SEQLEN], head[400];
	int m, i, j, k, v, fl = 0;
	int n_backr[MOTLEN];
	double freq[16];
	FILE  *in;

	int len_max = 0;
	for (k = 0; k < nseq; k++)
	{
		if (len[k] > len_max)len_max = len[k];
	}
	int *seq_back;
	seq_back = new int[len_max - 1];
	if (seq_back == NULL) { puts("Out of memory..."); exit(1); }
	for (j = 1; j <= reg_max; j++)
	{
		int j1 = j - 1;
		n_backr[j1] = 0;
		for (k = 0; k < nseq; k++)n_backr[j1] += (len[k] - j);
		n_backr[j1] *= 2;
	}
	if ((in = fopen(file, "rt")) == NULL)
	{
		printf("ReadSeq! Input file %s can't be opened!\n", file);
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
			int lenx = strlen(d);
			int check = CheckStr(file, d, nn, 0);
			nn++;
			if (lenx >= olen && check == 1)
			{
				TransStr(d);
				d[len[n]] = '\0';
				int cod[2];
				for (m = 0; m < 2; m++)
				{
					if (m == 1)ComplStr(d);
					cod[0] = IdeLet(d[0]);
					for (k = 0; k < len[n] - 1; k++)
					{
						cod[1] = IdeLet(d[k + 1]);
						int code = 4 * cod[0] + cod[1];
						seq_back[k] = code;
						cod[0] = cod[1];
					}
					for (j = 1; j <= reg_max; j++)
					{
						int j1 = j - 1;
						for (k = 0; k < len[n] - j; k++)
						{
							for (i = 0; i < 16; i++)freq[i] = 0;
							for (v = 0; v < j; v++)freq[seq_back[k + v]]++;
							for (i = 0; i < 16; i++)
							{
								if (freq[i] > 0)
								{
									freq[i] /= j;
									dav[j1][i] += freq[i];
									dcv[j1][i] += freq[i] * freq[i];
								}
							}
						}
					}
				}
				n++;
			}
			else
			{
				if (lenx < olen)printf("peak %d (Len %d) ignored\n", n + 1, lenx);
				if (check == -1)printf("Unusual symbol, peak %d ignored\n%s\n", n + 1, d);
			}
			if (fl == -1)
			{
				for (k = 0; k < 16; k++)
				{
					for (j = 0; j < reg_max; j++)
					{
						dav[j][k] /= n_backr[j];
						dcv[j][k] /= n_backr[j];
						dcv[j][k] -= dav[j][k] * dav[j][k];
					}
				}
				fclose(in);
				delete[] seq_back;
				break;
			}
		}
		if (*l == symbol)
		{
			memset(head, 0, sizeof(head));
			DelHole(l);
			strcpy(head, l);
			fl = 0; continue;
		}
		if (fl == 0)
		{
			memset(d, 0, sizeof(d));
			DelHole(l);
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
		DelHole(l);
		strcat(d, l);
	}
}
int main(int argc, char *argv[])
{
	int *len, nseq, *lenb, nseqb, i, j, k, m;
	char file_for[500], file_back[500], path_fasta[500], pfile_for[500], pfile_back[500];
	char filef[500], file_out[500];
	int ***seq_real;
	char ***peak_real;
	double **dav;//dinucl.content background
	double **dcv;//self covariations for regions LPD
	double *best_sco;
	double **frp;//LPD frequencies

	if (argc != 7)
	{
		puts("Sintax: 1path_both_fasta 2file_forground 3file_background 4int max_LPD_length 5int motif_len 6int size");//  5<pop_size>
		exit(1);
	}
	strcpy(path_fasta, argv[1]);
	strcpy(file_for, argv[2]);
	strcpy(file_back, argv[3]);
	strcpy(pfile_for, path_fasta);
	strcpy(pfile_back, path_fasta);
	strcat(pfile_for, file_for);
	strcat(pfile_back, file_back);
	int reg_max = atoi(argv[4]);// max LPD length
	int olen = atoi(argv[5]);// dlina motiva
	int size = atoi(argv[6]);
	int olen1 = olen - 1;
	srand((unsigned)time(NULL));
	if (size>POPSIZE)
	{
		printf("Maximal population size too large %d\n", size);
		exit(1);
	}
	dcv = new double *[reg_max];
	if (dcv == NULL) return -1;
	for (i = 0; i < reg_max; i++)
	{
		dcv[i] = new double[16];
		if (dcv[i] == NULL) return -1;
	}
	dav = new double *[reg_max];
	if (dav == NULL) return -1;
	for (i = 0; i < reg_max; i++)
	{
		dav[i] = new double[16];
		if (dav[i] == NULL) return -1;
	}
	town_ext pop_ext;
	//foreground
	EvalSeq(pfile_for, nseq, olen);
	len = new int[nseq];
	if (len == NULL) { puts("Out of memory..."); exit(1); }
	EvalLen(pfile_for, len, olen);
	//background
	EvalSeq(pfile_back, nseqb, olen);
	lenb = new int[nseqb];
	if (lenb == NULL) { puts("Out of memory..."); exit(1); }
	EvalLen(pfile_back, lenb, olen);
	best_sco = new double[nseq];
	if (best_sco == NULL){ puts("Out of memory..."); exit(1); }
	seq_real = new int**[2];
	if (seq_real == NULL){ puts("Out of memory..."); exit(1); }
	for (i = 0; i<2; i++)
	{
		seq_real[i] = new int*[nseq];
		for (j = 0; j < nseq; j++)
		{
			seq_real[i][j] = new int[len[j] - 1];
			if (seq_real[i][j] == NULL){ puts("Out of memory..."); exit(1); }
		}
	}
	peak_real = new char**[2];
	if (peak_real == NULL){ puts("Out of memory..."); exit(1); }
	for (i = 0; i<2; i++)
	{
		peak_real[i] = new char*[nseq];
		for (j = 0; j < nseq; j++)
		{
			peak_real[i][j] = new char[len[j] + 1];
			if (peak_real[i][j] == NULL){ puts("Out of memory..."); exit(1); }
		}
	}
	for (j = 0; j < reg_max; j++)
	{
		for (i = 0; i < 16; i++)dcv[j][i] = dav[j][i] = 0;
	}
	ReadSeq(pfile_for, nseq, len, seq_real, peak_real, olen);
	ReadSeqBack(pfile_back, nseqb, lenb, olen, reg_max, dav, dcv);
	for (j = 0; j < reg_max; j++)
	{
		printf("AvR%d\t", j + 1);
		for (i = 0; i < 16; i++)printf("%.4f\t", dav[j][i]);
		printf("\n");
	}
	for (j = 0; j < reg_max; j++)
	{
		printf("CovR%d\t", j + 1);
		for (i = 0; i < 16; i++)printf("%.4f\t", dcv[j][i]);
		printf("\n");
	}
	{		
		frp = new double *[nseq];
		if (frp == NULL) { puts("Out of memory..."); exit(1); }
		for (k = 0; k < nseq; k++)
		{
			frp[k] = new double[size];
			if (frp[k] == NULL) { puts("Out of memory..."); exit(1); }
		}
		for (i = 0; i < nseq; i++)
		{
			for (j = 0; j < size; j++)frp[i][j] = 0;
		}
	}
	{
		char word[] = "acgt";
		GetWords(2, 0, 16, word);
	}
	strcpy(filef, file_for);
	strcat(filef, "_");
	strcpy(file_out, file_for);
	strcat(file_out, "_andy");
	for (i = 0; i<MEGE; i++)
	{
		pop[i].mem_in(nseq);
	}
	for (i = 0; i<2; i++)
	{
		det2[i].mem_in(nseq);
	}
	det1.mem_in(nseq);

	{
		int big_exit1 = 1;// local exit (separ +-) global exit (separation do not exceeded the previous run)
		double fit_prev, fit_after_mut;
		int cnt_count = 0;
		//Test(peak_real[0],len,0,2);		
		printf("\n%s\tTrain\t", file_for);
		printf("Ndi %d\tDeg %d\tEli %d BE1 %d\n", size, MEGE, ELIT, big_exit1);
		//initiation		
		int gen = 0;
		int restart = 0;
		int rec_first_only = 0;
		int pair_all = MEGE*(MEGE - 1);
		int pair_d[MEGE*(MEGE - 1)][2];
		int pair_take[MEGE*(MEGE - 1)];
		int success_r, success_m;
		int success_r1[NREC];
		//PARAMETERS SETTING																
		int stop_oi[MEGE], stop_li[MEGE], stop_pi[MEGE];
		for (i = 0; i < MEGE; i++)stop_oi[i] = stop_li[i] = stop_pi[i] = 0;
		int mege_h;
		do
		{
			int success_o, success_l, success_p;
			if (big_exit1 == 1)
			{
				//if(isize==0)
				{
					for (i = 0; i < MEGE; i++)
					{
						pop[i].fit = 0;
						town ini;
						ini.mem_in(nseq);
						for (j = 0; j < MEGE; j++)
						{
							int gom = 0;
							do
							{
								int m1;
								do
								{
									gom = 0;
									ini.init_rand(nseq, len, olen, size, reg_max);
									if (ini.check(0, reg_max) == -1)
									{
										ini.check(0, reg_max);
										printf("Population error!\n");
										ini.print_all(reg_max, nseq);
										gom = -1;//exit(1);
									}
									if (gom == 0)
									{
										for (m = 0; m < i; m++)
										{
											gom = GomTown(ini, pop[m], nseq);
											if (gom == -1)
											{
												m1 = m;
												break;
											}
										}
									}
								} while (gom == -1);
								EvalMahFIT(&ini, nseq, seq_real, dav, dcv,frp);
								ini.print_sta(reg_max);
								if (gom == 0)
								{
									ini.get_copy(&pop[i], nseq, reg_max);
								}
							} while (ini.fit == 0 || gom == -1);
							if (ini.fit > 0)break;
						}
						//qsort((void*)(&ini[0]),n_rand, sizeof(ini[0]), compare_pop);					
						ini.mem_out();
					}
				}
				big_exit1 = 0;
				qsort((void*)(&pop[0]), MEGE, sizeof(pop[0]), compare_pop);
			}
			/*for (i = 0; i<MEGE; i++)
			{
			//printf("After %d\n",i+1);
			//	pop[i].print_all(reg_max,nseq);
			if (pop[i].check(0, reg_max) == -1)
			{
			pop[i].print_all(reg_max, nseq);
			printf("Population error!\n");
			exit(1);
			}
			}*/
			fit_prev = pop[0].fit;
			//pop[0].print_all(reg_max,nseq);													
			success_o = success_l = success_p = success_m = 0;
			double ratio_thr = 0.004, ratio_thr_r = ratio_thr;
			int step;
			int step_max, step_max_tot = 0;
			if (restart == 0)
			{
				step_max = 1000;
				step = 200;
				mege_h = MEGE;
			}
			else
			{
				//step_max = 10000 * (int)(nseq / 20);
				step_max = 100000;//10000 * (int)(nseq / 20);
				step = 10000;
				mege_h = ELIT;
			}
			int mege_h1 = mege_h - 1;
			int step2 = 2 * step;
			//mutations						
			int n_mut_tot = 0, success_m_tot = 0, mdo = 1;
			int asuccess[NMUT], atry[NMUT];
			for (k = 0; k < NMUT; k++)asuccess[k] = atry[k] = 0;
			int success_mi[MEGE], try_mi[MEGE];
			int success_ri[MEGE];
			double exp_rec_rate[MEGE];
			for (i = 0; i < mege_h; i++)exp_rec_rate[i] = 0;
			for (i = 0; i < mege_h; i++)success_mi[i] = try_mi[i] = 0;
			do
			{
				step_max_tot += step_max;
				if (restart == 0)printf("Mut cycle %d\n", step_max_tot / step_max);
				for (i = 0; i < mege_h; i++)
				{
					if (stop_pi[i] == 0 || (stop_li[i] == 0 || stop_oi[i] == 0))
					{
						//	printf("Mut start %d\n",i+1);			
						pop[i].get_copy(&det1, nseq,reg_max);
						int success_o_local = 0, success_l_local = 0, success_p_local = 0;
						int success_m_local = 0;
						double ratio[NMUT];
						if (stop_oi[i] == 0)ratio[0] = 1;
						else ratio[0] = 0;
						if (stop_li[i] == 0)ratio[1] = 1;
						else ratio[1] = 0;
						if (stop_pi[i] == 0)ratio[2] = 1;
						else ratio[2] = 0;
						int step_success[NMUT];
						int step_try[NMUT];
						for (k = 0; k < NMUT; k++)step_success[k] = step_try[k] = 0;
						int n_mut_here = 0;
						int n_mut_min_cyc = step;
						int m_iter = 0;
						while (ratio[2] > ratio_thr || (ratio[0] > ratio_thr || ratio[1] > ratio_thr))
						{
							//mm++;
							//printf("I=%d Mut %d\t%f %f %f\t",i+1,mm,ratio[0],ratio[1],ratio[2]);							
							int sm;
							//if (ratio[0] > ratio_thr && (ratio[1] > ratio_thr && ratio[2] > ratio_thr))
							//if (stop_li[i] == 1 && stop_oi[i] == 1)sm = 2;
							//	else
							{
								int rs = step + step_success[2];
								if (ratio[1] > ratio_thr && stop_li[i] == 0)rs += step + step_success[1];
								if (ratio[0] > ratio_thr && stop_oi[i] == 0)rs += step + step_success[0];
								int rs1 = rand() % rs;
								if (rs1 < step + step_success[2])sm = 2;
								else
								{
									if (rs1 < step2 + step_success[2] + step_success[1])sm = 1;
									else sm = 0;
								}
								//sm = rand() % 3;
							}
							int muto, muto1;
							int npeak = -1, nori = -1, npos = -1;
							if (sm == 0)muto = MutOlig0(&det1);
							else
							{
								if (sm == 1)muto = MutCry0(&det1, olen, reg_max);
								else
								{
									muto = MutRegShift(&det1, nseq, len, olen, npeak, nori, npos);
								}
							}
							int gom = 0;													
							if (muto != -1)
							{
								muto1 = det1.order(muto);
								if (muto1 == -1)
								{
									puts("MUTO Error");
									gom = -1;
								}
								if (gom == 0)
								{
									for (m = 0; m < mege_h; m++)
									{
										if (m == i)continue;
										gom = GomTown(det1, pop[m], nseq);
										if (gom == -1) { break; }
									}
								}
							}									
							//if (det1.check(0, reg_max) == -1){ printf("Population error!\n"); exit(1); }
							//det.print_all();										
							if (gom != -1)
							{
								double dd = 0;
								n_mut_here++;
								EvalMahFIT(&det1, nseq, seq_real, dav, dcv,frp);
								dd = det1.fit / pop[i].fit;
								if (dd > 1)
								{
									step_success[sm]++;
									if (sm == 0)success_o_local++;
									else
									{
										if (sm == 1)success_l_local++;
										else success_p_local++;
									}
									success_m_local++;
									success_mi[i]++;
									det1.get_copy(&pop[i], nseq, reg_max);
								}
								else
								{
									pop[i].get_copy(&det1, nseq, reg_max);
								}
								step_try[sm]++;
								try_mi[i]++;
								if (n_mut_here % step == 0)
								{
									m_iter++;
									for (k = 0; k < 3; k++)
									{
										if (step_success[k] != 0 && step_try[k] != 0)
										{
											ratio[k] = (double)step_success[k] / step_try[k];
										}
										else ratio[k] = 0;
									}
									//		for (k = 0; k < 3; k++)printf("%f ", ratio[k]); printf("\t");
									//	for (k = 0; k < 3; k++)printf("%d ", step_success[k]); printf("\t");
									//for (k = 0; k < 3; k++)printf("%d ", step_try[k]); printf("\t");
									//	double ratio_thr_m1 = (double)(step_success[0] + step_success[1] + step_success[2]) / (step_try[0] + step_try[1] + step_try[2]);									
									int step_success_tot = step_success[0] + step_success[1] + step_success[2];
									if (step_success_tot < n_mut_min_cyc)n_mut_min_cyc = step_success_tot;									
									if (restart == 1)printf("Step%d M%d %d Try %d Min %d Mah %f Ave %f Std %f Fit %f Ratio %f\n", m_iter, i + 1, step_success[2], step_try[2], n_mut_min_cyc, pop[i].mah, pop[i].ave, pop[i].std, pop[i].fit, ratio[2]);
									if (ratio[0] <= ratio_thr)stop_oi[i] = 1;
									if (ratio[1] <= ratio_thr)stop_li[i] = 1;
									if (ratio[2] <= ratio_thr)
									{
										stop_pi[i] = 1;
										//step_max = m_iter*step;
										//step_max = n_mut_tot;
									}
									for (k = 0; k < 3; k++)atry[k] += step_try[k];
									for (k = 0; k < 3; k++)asuccess[k] += step_success[k];
									for (k = 0; k < 3; k++)step_try[k] = step_success[k] = 0;
									if (try_mi[i] != 0)exp_rec_rate[i] = (double)success_mi[i] / try_mi[i];
									else exp_rec_rate[i] = 0;
									success_mi[i] = try_mi[i] = 0;
									if (n_mut_here >= step_max)break;
									//if (restart==0 && n_mut_here >= step_max)break;									
									//if (restart == 1 && stop_pi[i] == 1)break;
								}
								//printf("\n");
							}
						}
						success_o += success_o_local;
						success_l += success_l_local;
						success_p += success_p_local;
						success_m_tot += success_m_local;
						printf("M %d %d,%d,%d = %d Try %d Min %d M %f A %f S %f F %f", i + 1, success_o_local, success_l_local, success_p_local, success_m_local, n_mut_here, n_mut_min_cyc, pop[i].mah, pop[i].ave, pop[i].std, pop[i].fit);
						//if (restart == 0)
						printf("\tOLP %d%d%d", stop_oi[i], stop_li[i], stop_pi[i]);
						printf("\n");
						n_mut_tot += n_mut_here;
						if (success_m_local != 0)success_m++;
					}
				}
				if (restart == 1)mdo = 0;
				else
				{
					mdo = 0;
					for (i = 0; i < mege_h; i++)
					{
						if (stop_li[i] == 0 || stop_oi[i] == 0)
						{
							mdo = 1;
							break;
						}
					}
				}
			} while (mdo == 1);
			qsort((void*)(&pop[0]), mege_h, sizeof(pop[0]), compare_pop);
			if (stop_pi[0] == 1)
			{
				rec_first_only = 1;
				//		printf("Rec First only!\n");
			}
			if (restart == 0)success_m /= (step_max_tot / step_max);
			if (restart == 0)ratio_thr_r = (double)(asuccess[0] + asuccess[1]) / (atry[0] + atry[1]);
			else
			{
				//	if (stop_pi[0] == 1)
				{
					int recount = 0;
					double resum = 0;
					for (i = 0; i < mege_h1; i++)
					{
						double e1 = Max(ratio_thr, exp_rec_rate[i]);
						double e2 = Max(ratio_thr, exp_rec_rate[i + 1]);
						resum += success_ri[i] * (e1 + e2) / 2;
						recount += success_ri[i];
					}
					if (recount != 0)ratio_thr_r = resum / recount;
					else ratio_thr_r = ratio_thr;
				}
				//	else ratio_thr_r = exp_rec_rate[mege_h1];
			}
			//recombinations											
			printf("Total mut %d\n", n_mut_tot);
			int loc_rec, loc_rec_tot = 0;
			fit_after_mut = pop[0].fit;
			for (m = 0; m<NREC; m++)success_r1[m] = 0;
			pair_all = 0;
			int jmax;
			if (rec_first_only == 0)
			{
				if (gen> 0)jmax = mege_h / 3 - 1;
				else jmax = mege_h / 2 - 1;
			}
			else
			{
				if (gen > 1)jmax = mege_h / 3 - 1;
				else jmax = mege_h / 2 - 1;
			}
			{
				int multw = 1;
				for (j = jmax; j >= 0; j--)
				{
					for (k = j + 1; k < mege_h; k++)
					{
						for (m = 0; m < multw; m++)
						{
							pair_d[pair_all][0] = k;
							pair_d[pair_all][1] = j;
							//							printf("%d, %d\t", k, j);
							pair_all++;
							pair_d[pair_all][0] = j;
							pair_d[pair_all][1] = k;
							//						printf("%d, %d\n", j, k);
							pair_all++;
						}
					}
					if (gen > 1)multw++;
				}
			}
			for (k = 0; k<pair_all; k++)pair_take[k] = k;
			int n_rec_cycle, sr_step;
			{
				int max_rec = Max(step_max, n_mut_tot);
				n_rec_cycle = max_rec / pair_all;
				if (gen == 0)max_rec = step_max_tot;
				else max_rec = step;
				sr_step = max_rec / pair_all;
			}
			if (sr_step < 1)sr_step = 1;
			if (n_rec_cycle<1)n_rec_cycle = 1;
			int sr;
			double ratio_r = 1;
			int step_rtry[NREC], step_rsuccess[NREC];
			for (k = 0; k < NREC; k++)step_rtry[k] = step_rsuccess[k] = 0;
			int success_r_cycle = 0;
			printf("Rec %d cycles of %d tries\tStep %d\tRatioThr %.4f RecOsob 0..%d\n", n_rec_cycle, pair_all, step_max_tot, ratio_thr_r, jmax);
			loc_rec = 0;
			int elit_rec = ELIT;
			double fit_rec_prev = pop[elit_rec].fit;
			for (i = 0; i < mege_h; i++)success_ri[i] = 0;
			for (sr = 1; sr <= n_rec_cycle; sr++)
			{
				BigMixI(pair_take, pair_all);
				for (k = 0; k<pair_all; k++)
				{
					int kk[2];
					for (m = 0; m<2; m++)kk[m] = pair_d[pair_take[k]][m];
					//if (rec_first_only == 1 && kk[1] != 0)continue;
					for (m = 0; m<2; m++)
					{
						pop[kk[m]].get_copy(&det2[m], nseq, reg_max);
					}
					double fit_parent_max;
					fit_parent_max = Max(det2[0].fit, det2[1].fit);
					int nsn[2], nsi[2], num[2];
					int sumw = 0, r_cy;
					int kk01 = pair_take[k];
					if (gen == 0)r_cy = rand() % 5;
					else
					{
						if (restart == 0)r_cy = rand() % 5;
						else r_cy = 4;
					}
					switch (r_cy)
					{
					case(0) :
					{
						if (Reco2_Original(&det2[0], &det2[1], nsi, num) == -1)continue;
						else break;
					}
					case(1) :
					{
						if (Reco2_Economic(&det2[0], &det2[1], nsi, num) == -1)continue;
						else break;
					}
					case(2) :
					{
						if (Reco2_One_dinucleotide_local(&det2[0], &det2[1], olen, reg_max) == -1)continue;
						else break;
					}
					case(3) :
					{
						if (Reco2_One_dinucleotide_full(&det2[0], &det2[1]) == -1)continue;
						else break;
					}
					case(4) :
					{
						if (Reco2Peak(&det2[0], &det2[1], nseq) == -1)continue;
						else break;
					}
					}
					//printf("out Rec %d\n",r_cy);
					if (r_cy <= 1)
					{
						for (m = 0; m<2; m++)
						{
							//	det[m].tot[nsi[m]].print_all();
							nsn[m] = det2[m].order(nsi[m]);
						}
					}
					/*for (m = 0; m<2; m++)
					{
					if (det2[m].check(0, reg_max) == -1)
					{
					det2[m].check(0, reg_max);
					printf("Population error!\n");
					exit(1);
					}
					}*/
					double dd_r[2] = { 0, 0 };
					int gom[2] = { 0, 0 };
					{
						for (m = 0; m<2; m++)
						{
							for (int t = 0; t< mege_h; t++)
							{
								if (t == pair_d[k][m])continue;
								gom[m] = GomTown(det2[m], pop[t], nseq);
								if (gom[m] == -1)
								{
									/*	if (pop[t].fit>fit_parent_max)
									{
									fit_parent_max = pop[t].fit;
									pair_d[k][m] = t;
									kk[m] = pair_d[pair_take[k]][m];
									kk01 = t;
									}
									else gom[m] = 1;*/
									break;
								}
							}
							if (gom[m] == -1)break;
						}
					}
					if (gom[0] != -1 && gom[1] != -1)
					{
						for (m = 0; m < 2; m++)
						{
							det2[m].fit = EvalMahFIT(&det2[m], nseq, seq_real, dav, dcv,frp);
						}
						/*
						if(r_cy==4)
						{
						printf("Po\n");
						det[0].print_all(reg_max,nseq);
						det[1].print_all(reg_max,nseq);
						}*/
						step_rtry[r_cy]++;
						double fit_det_max = Max(det2[0].fit, det2[1].fit);
						if (fit_det_max > fit_parent_max)
						{
							int kk_min = Min(kk[0], kk[1]);
							success_ri[kk_min]++;
							//success_ri[kk[0]]++;
							//success_ri[kk[1]]++;
							success_r1[r_cy]++;
							step_rsuccess[r_cy]++;
							if (restart == 0)
							{
								if (r_cy <= 3 && fit_det_max > fit_rec_prev)loc_rec++;
							}
							success_r_cycle++;
							for (m = 0; m < 2; m++)
							{
								det2[m].get_copy(&pop[kk[m]], nseq, reg_max);
							}
						}
					}
				}
				success_r = 0;
				for (i = 0; i< mege_h; i++)if (success_ri[i]>0)success_r++;
				if (sr%sr_step == 0)
				{
					if (restart == 1)ratio_r = (double)success_r_cycle / pair_all / sr_step;
					else
					{
						int stepr = (step_rtry[0] + step_rtry[1] + step_rtry[2]);
						if (stepr > 0)ratio_r = (double)(step_rsuccess[0] + step_rsuccess[1] + step_rsuccess[2]) / stepr;
						else ratio_r = 0;
					}
					int mbe = 0;
					for (m = 1; m < mege_h; m++)
					{
						if (pop[m].fit > pop[0].fit)
						{
							mbe = m;
						}
					}					
					//					printf("Rcycle %d: %d %d,%d,%d,%d,%d %d %f %f\n", sr, success_r, success_r1[0], success_r1[1], success_r1[2], success_r1[3], success_r1[4], success_r_cycle, ratio_r, pop[0].fit);
					if (mbe != 0)qsort((void*)(&pop[0]), mege_h, sizeof(pop[0]), compare_pop);
					fit_rec_prev = pop[elit_rec].fit;
					loc_rec_tot += loc_rec;
					printf("Rec %d: %d %d,%d,%d,%d,%d %d %f M %f A %f S %f F %f ", sr*pair_all, success_r, success_r1[0], success_r1[1], success_r1[2], success_r1[3], success_r1[4], success_r_cycle, ratio_r, pop[mbe].mah, pop[mbe].ave, pop[mbe].std, pop[mbe].fit);
					if (restart == 0)printf("L%d", loc_rec_tot);
					printf("\n");
					if (gen > 0 && loc_rec == 0)break;
					if (ratio_r < ratio_thr_r)
					{
						if (restart == 0)
						{
							if (loc_rec == 0)break;
						}
						else break;
					}
					loc_rec = 0;
					for (m = 0; m < NREC; m++)step_rtry[m] = step_rsuccess[m] = 0;
					success_r_cycle = 0;
				}
			}
			qsort((void*)(&pop[0]), mege_h, sizeof(pop[0]), compare_pop);
			gen++;
			double change_level = pop[0].fit / fit_prev - 1;
			double change_level_rec = pop[0].fit / fit_after_mut - 1;
			double change_level_mut = fit_after_mut / fit_prev - 1;
			if (restart == 0)
			{
				restart = 1;
				for (i = 1; i < ELIT; i++)
				{
					pop[0].get_copy(&pop[i], nseq, reg_max);
				}
				for (i = 0; i < ELIT; i++)
				{
					pop[i].init_rand_part(nseq, len, olen, 20);
					EvalMahFIT(&pop[i], nseq, seq_real, dav, dcv,frp);
				}
				qsort((void*)(&pop[0]), ELIT, sizeof(pop[0]), compare_pop);
			}
			fit_prev = pop[0].fit;
			//if (change_level<GA_EXIT)gen1++;
			printf("Gen %d Fit %.5f Rat %.5f RatM %.5f RatR %.5f ", gen, pop[0].fit, change_level, change_level_mut, change_level_rec);
			//	if(restart==1)printf("W %d WP %d AvS %.1f AvP %.1f", wei_max, wei_max_pos, wei_peak_av, wei_pos_av);
			printf("M %d %d,%d,%d R %d ", success_m, success_o, success_l, success_p, success_r1[4]);
			//	if (restart == 0)printf("LM %d LR %d ", loc_mut_tot, loc_rec_tot);
			{
				//	printf("R %d %d,%d,%d,%d,%d ", success_r, success_r1[0], success_r1[1], success_r1[2], success_r1[3], success_r1[4]);
				//	printf("%d,%d,%d,%d,%d ", recw[0][0], recw[1][0], recw[2][0], recw[3][0], recw[4][0]);
				//	printf("%d,%d,%d,%d,%d ", recw[0][pair_all1], recw[1][pair_all1], recw[2][pair_all1], recw[3][pair_all1], recw[4][pair_all1]);				
			}
			{
				int sumr = 0;
				for (i = 0; i < mege_h; i++)sumr += success_ri[i];
				for (i = 0; i < mege_h; i++)printf("%.2f ", 100 * (double)success_ri[i] / sumr);
			}
			{
				time_t tnow;
				time(&tnow);
				printf("%s", ctime(&tnow));
			}
			fit_prev = pop[0].fit;
			if (restart == 0)for (i = 0; i<2; i++)pop[i].print_all(reg_max, nseq);
			if (stop_pi[0] == 1 && gen > 1)
			{
				big_exit1 = 1;
				printf("Go out ");
			}
		} while (big_exit1 == 0);
		EvalMahFITTrain(&pop[0], nseq, seq_real, file_for, olen, len, peak_real, &pop_ext, best_sco, dav, dcv);
		printf("Go out big cycle ");
		big_exit1 = 1;
	}
	char name[500];
	for (i = 0;; i++)
	{
		if (file_for[i] == '.'){ name[i] = '\0'; break; }
		if (file_for[i] == '\n'){ name[i] = '\0'; break; }
		if (file_for[i] == '\0'){ name[i] = '\0'; break; }
		name[i] = file_for[i];
	}
	{
		char ext2best[20], extmat[20];
		char ext2best0[] = "_end2";
		char extmat0[] = "_mat";
		strcpy(ext2best, ext2best0);
		strcpy(extmat, extmat0);
		//	strcat(ext2best, argv[2]);
		//	strcat(extmat, argv[2]);
		pop[0].fprint_allfi(file_for, ext2best, olen, pop_ext.c0, pop_ext.buf);
		pop[0].fprint_allfi_mat(file_for, extmat, name, olen, pop_ext.c0, pop_ext.buf);
	}
	char file_train_seq[500];
	int file_fasta_len = strlen(file_for);
	k = 0;
	for (j = 0; j<file_fasta_len; j++)
	{
		char cc = file_for[j];
		if (cc == '.' || cc == '\0')
		{
			file_train_seq[k++] = '\0';
			break;
		}
		file_train_seq[k++] = cc;
	}
	strcat(file_train_seq, "_sga_train.seq");
	pop[0].fprint_seq(file_train_seq, olen, nseq, peak_real, best_sco);
	for (i = 0; i<MEGE; i++)pop[i].mem_out();
	for (i = 0; i<2; i++)det2[i].mem_out();
	det1.mem_out();
	for (k = 0; k<2; k++)
	{
		for (i = 0; i<nseq; i++)
		{
			delete[] seq_real[k][i];
		}
		delete[] seq_real[k];
	}
	delete[] seq_real;
	for (k = 0; k<2; k++)
	{
		for (i = 0; i<nseq; i++)
		{
			delete[] peak_real[k][i];
		}
		delete[] peak_real[k];
	}
	delete[] peak_real;
	delete[] best_sco;
	delete[] len;
	delete[] lenb;
	for (i = 0; i < reg_max; i++)
	{
		delete[] dcv[i];
	}
	delete[] dcv;
	for (i = 0; i < reg_max; i++)
	{
		delete[] dav[i];
	}
	delete[] dav;
	for (k = 0; k < nseq; k++)
	{
		delete[] frp[k];
	}
	return 0;
}
