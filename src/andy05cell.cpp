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
#define LPDLEN 6 //max LPD length
#define MOTLEN 20 //max motif len
#define MEGE 100//population size 1st stage
#define NMUT 3
#define NREC 6
#define POPSIZE 120
//#define GA_EXIT 0.01

double  uw[POPSIZE][POPSIZE];
struct ss {
	int num;
	char oli[3];
}s[16];

int compare_qq(const void* X1, const void* X2)
{
	double X = (*(double*)X1 - *(double*)X2);
	if (X > 0)return 1;
	if (X < 0)return -1;
	return 0;
}
int compare_qq2(const void* X1, const void* X2)
{
	double X = (*(double*)X1 - *(double*)X2);
	if (X > 0)return -1;
	if (X < 0)return 1;
	return 0;
}
struct qbs {
	double q;
	int n;
};
int compare_qbs(const void* X1, const void* X2)
{
	struct qbs* S1 = (struct qbs*)X1;
	struct qbs* S2 = (struct qbs*)X2;
	if (S1->q - S2->q > 0)return 1;
	if (S1->q - S2->q < 0)return -1;
	return 0;
}
struct uno {
	int sta;
	int end;
	int num;
	void get_copy(uno* a);
	void print_all(FILE* outlog);
};
void uno::get_copy(uno* a)
{
	a->num = num;
	a->sta = sta;
	a->end = end;
};
void uno::print_all(FILE* outlog)
{
	fprintf(outlog, "[%d;%d]%s\t", sta, end, s[num].oli);
}
struct town {
	uno tot[POPSIZE];
	int deg[16];
	int size;
	int* pos;// pozicii na4al okon
	int* ori;// DNA strand 0,1
	double mah;
	double fpr;
	double inf;
	double fit;
	int sm;
	int tm;
	int sr;
	int tr;
	int odg[LPDLEN + 2];
	void get_copy(town* a, int nseq, int reg_max);
	//void init_rand(int nseq, int *len, int oln, int rsize, int reg_max);
	void init_rand_hoxa(int nseq, int oln, int rsize, int reg_max, int* len_octa, int* len, int** octa_prowb, int** octa_prows);
	void init_rand_part_hoxa(int nseq, int nind, int olen, int* len_octa, int* len, int** octa_prowb, int** octa_prows);
	int init_add(uno last);
	int init_add_split(void);
	void init_zero(int olen);
	void Mix(int* a, int* b);
	void BigMix(int* d, int ln);
	void swap(int n1, int n2);
	int order(int n);
	void print_all(int reg_max, int nseq, FILE* outlog);
	void print_sta(int reg_max, FILE* outlog);
	void sort_all(void);
	int sum(int j);
	void fprint_all(char* file, char* add);
//	void fprint_allfi(char* file, char* add, int len, double c0, double* buf, int reg_max);
	void fprint_allfi_mat(char* file, char* add, char* name, int len, double sga_min, double sga_raz, double* buf);	
	int check(int min, int max, FILE* out);
	int mem_in(int nseq);
	void mem_out(void);
} pop[MEGE], det1, det2[2];
int compare_tot(const void* X1, const void* X2)
{
	struct uno* S1 = (struct uno*)X1;
	struct uno* S2 = (struct uno*)X2;
	if (S1->num - S2->num > 0)return 1;
	if (S1->num - S2->num < 0)return -1;
	if (S1->sta - S2->sta > 0)return 1;
	if (S1->sta - S2->sta < 0)return -1;
	return 0;
}
int compare_pop(const void* X1, const void* X2)
{
	struct town* S1 = (struct town*)X1;
	struct town* S2 = (struct town*)X2;
	if (S1->fit - S2->fit > 0)return -1;
	if (S1->fit - S2->fit < 0)return 1;
	return 0;
}
void town::init_zero(int reg_max)
{
	int k;
	for (k = 0; k < 16; k++)deg[k] = 0;
	for (k = 0; k < reg_max; k++)odg[k] = 0;
	odg[reg_max] = -1;
	size = 0;
	fit = mah = 0;
	fpr = inf = 1;
}
int town::init_add(uno last)
{
	int k;
	int dtype = last.num;
	int sta = last.sta;
	int end = last.end;
	for (k = 0; k < size; k++)
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
int town::check(int min, int max, FILE* out)
{
	int i, j;
	int ods[LPDLEN];
	int max1 = max - 1;
	for (i = 0; i < max1; i++)ods[i] = 0;
	for (i = 0; i < size; i++)
	{
		ods[tot[i].end - tot[i].sta]++;
	}
	for (i = 0; i < max1; i++)
	{
		if (ods[i] != odg[i])
		{
			fprintf(out, "odg(%d) %d actual %d", i, odg[i], ods[i]);
			return -1;
		}
	}
	int sum = 0;
	for (i = 0; i < max; i++)
	{
		if (odg[i] < 0) { fprintf(out, "odg %d = %d", i, odg[i]); return -1; }
		sum += odg[i];
	}
	if (sum != size)
	{
		fprintf(out, "wrong sum odg = %d, size = %d", sum, size);
		return -1;
	}
	for (i = 0; i < size; i++)
	{
		int len = tot[i].end - tot[i].sta;
		if (len < min || len >= max)
		{
			tot[i].print_all(out);
			fprintf(out, "\n");
			return -1;
		}
	}
	for (i = 0; i < size; i++)
	{
		for (j = i + 1; j < size; j++)
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
					fprintf(out, "1st %d\t", i + 1);
					tot[j].print_all(out);
					fprintf(out, "\n");
					fprintf(out, "2nd %d\t", j + 1);
					tot[i].print_all(out);
					fprintf(out, "\n");
					return -1;
				}
			}
		}
	}
	return 1;
}
void MixI(int* a, int* b)
{
	int buf = *a;
	*a = *b;
	*b = buf;
}
void BigMixI(int* d1, int len) // pereme6ivanie stroki
{
	int r;
	for (r = 0; r < len - 1; r++)
	{
		MixI(&d1[r], &d1[1 + r + (rand() % (len - 1 - r))]);
	}
}
void MixC(char* a, char* b)
{
	char buf = *a;
	*a = *b;
	*b = buf;
}
void BigMixC(char* d1, int len) // pereme6ivanie stroki
{
	int r;
	for (r = 0; r < len - 1; r++)
	{
		MixC(&d1[r], &d1[1 + r + (rand() % (len - 1 - r))]);
	}
}
void town::print_all(int reg_max, int nseq, FILE* outlog)
{
	int i;
	char strand[] = "+-";
	fprintf(outlog, "M %f FP %f F %f\t", mah, fpr, fit);
	for (i = 0; i < 16; i++)fprintf(outlog, "%d ", deg[i]); fprintf(outlog, "\tLEN ");
	for (i = 0; i < reg_max; i++)fprintf(outlog, "%d ", odg[i]); fprintf(outlog, "\t");
	int size1 = size - 1;
	fprintf(outlog, " ");
	for (i = 0; i < size; i++)
	{
		tot[i].print_all(outlog);
		if (i == size1)fprintf(outlog, "\n");
		// else if(tot[i + 1].num != tot[i].num)fprintf(outlog,"\n");
		//	 fprintf(outlog," ");
	}
	// fprintf(outlog,"\n");
	/*
	for(i=0;i<nseq;i++)
	{
	fprintf(outlog,"%2d",pos[i]);
	fprintf(outlog,"%c",strand[ori[i]]);
	fprintf(outlog," ");
	}
	fprintf(outlog,"\n");*/
}
void town::print_sta(int reg_max, FILE* outlog)
{
	int i;
	fprintf(outlog, "M %f FP %f F %f\t", mah, fpr, fit);
	for (i = 0; i < 16; i++)fprintf(outlog, "%2d", deg[i]); fprintf(outlog, "\t");
	for (i = 0; i < reg_max; i++)fprintf(outlog, "%2d ", odg[i]); fprintf(outlog, "\n");
}
void town::get_copy(town* a, int nseq, int reg_max)
{
	int i;
	a->size = size;
	a->fit = fit;
	a->mah = mah;
	a->fpr = fpr;
	a->inf = inf;
	for (i = 0; i < size; i++)tot[i].get_copy(&a->tot[i]);
	for (i = 0; i < 16; i++)a->deg[i] = deg[i];
	i = 0;
	for (i = 0; i < reg_max; i++)
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
		if (tot[n].num < tot[n - 1].num)di = -1;
		else if (tot[n].num > tot[n + 1].num)di = 1;
		else
		{
			if (tot[n].num == tot[n - 1].num && tot[n].sta < tot[n - 1].sta)di = -1;
			if (tot[n].num == tot[n + 1].num && tot[n].sta > tot[n + 1].sta)di = 1;
		}
	}
	{
		if (n == 0)
		{
			if (tot[n].num > tot[n + 1].num)di = 1;
			else if (tot[n].num == tot[n + 1].num && tot[n].sta > tot[n + 1].sta)di = 1;
		}
		if (n == size - 1)
		{
			if (tot[n].num < tot[n - 1].num)di = -1;
			else if (tot[n].num == tot[n - 1].num && tot[n].sta < tot[n - 1].sta)di = -1;
		}
	}
	//printf("DIRECT %d\t",di);
	//if(di==0)return -1;
	if (di == -1)
	{
		for (i = 0; i <= n - 1; i++)
		{
			if (tot[n - i - 1].num > tot[n - i].num)
			{
				//			tot[n-i-1].print_all();
				swap(n - i - 1, n - i);
				ret--;
			}
			else if (tot[n - i - 1].num == tot[n - i].num && tot[n - i - 1].sta > tot[n - i].sta)
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
			if (tot[n + i + 1].num < tot[n + i].num)
			{
				//			tot[n+i+1].print_all();
				swap(n + i + 1, n + i);
				ret++;
			}
			else if (tot[n + i + 1].num == tot[n + i].num && tot[n + i + 1].sta < tot[n + i].sta)
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
void town::init_rand_hoxa(int nseq, int oln, int rsize, int reg_max, int* len_octa, int* len, int** octa_prowb, int** octa_prows)
{
	int i, j;

	fit = 0;
	fpr = 1;
	size = rsize;
	int oln1 = oln - 1;
	for (j = 0; j < reg_max; j++)odg[j] = 0;
	odg[reg_max] = -1;
	for (i = 0; i < nseq; i++)
	{
		ori[i] = rand() % 2;
		int inx = octa_prows[i][len_octa[i] - 1];
		int r2 = rand() % inx;
		for (j = 0; j < len_octa[i]; j++)
		{
			if (octa_prows[i][j] > r2)
			{
				pos[i] = octa_prowb[i][j];
				break;
			}
		}
	}
	for (i = 0; i < 16; i++)deg[i] = 0;
	i = 0;
	do
	{
		int r = rand() % 16;
		if (deg[r] == oln1)continue;
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
			for (j = deg[i]; j < oln1; j++)take_pos[j] = 0;
			//		printf("DO\t");
			//	for(j=0;j<oln1;j++)printf("%d\t",take_pos[j]);printf("\n");
			BigMixI(take_pos, oln1);
			//	printf("PO\t");
			//	for(j=0;j<oln1;j++)printf("%d\t",take_pos[j]);printf("\n");
			for (j = 0; j < oln1; j++)
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
/*void town::init_rand(int nseq, int *len, int oln, int rsize, int reg_max)
{
	int i, j;

	fit = 0;
	fpr = 1;
	size = rsize;
	int oln1 = oln - 1;
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
		if (deg[r] == oln1)continue;
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
			for (j = deg[i]; j < oln1; j++)take_pos[j] = 0;
			//		printf("DO\t");
			//	for(j=0;j<oln1;j++)printf("%d\t",take_pos[j]);printf("\n");
			BigMixI(take_pos, oln1);
			//	printf("PO\t");
			//	for(j=0;j<oln1;j++)printf("%d\t",take_pos[j]);printf("\n");
			for (j = 0; j < oln1; j++)
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
}*/
void town::init_rand_part_hoxa(int nseq, int nind, int olen, int* len_octa, int* len, int** octa_prowb, int** octa_prows)
{
	int i, k;

	fit = 0;
	fpr = 1;
	//int oln1 = olen - 1;
	for (i = 0; i < nseq; i++)
	{
		int r = rand() % nind;
		if (r == 0)
		{
			ori[i] = rand() % 2;
			int inx = octa_prows[i][len_octa[i] - 1];
			int r2 = rand() % inx;
			for (k = 0; k < len_octa[i]; k++)
			{
				if (octa_prows[i][k] > r2)
				{
					pos[i] = octa_prowb[i][k];
					break;
				}
			}
		}
	}
}
/*void town::init_rand_part(int nseq, int *len, int oln, int *xporti, int nind)
{
	int i, j;

	fit = 0;
	int oln1 = oln - 1;
	for (j = 0; j < nseq; j++)
	{
		i = xporti[j];
		int r = rand() % nind;
		if (r == 0)
		{
			int lenp = len[i] - oln1;
			pos[i] = rand() % lenp;
			ori[i] = rand() % 2;
		}
	}
}*/
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
	for (i = 0; i < j; i++)ret += deg[i];
	return ret;
}
void town::fprint_all(char* file, char* add)
{
	int i;
	FILE* out;
	char file_out[500];
	strcpy(file_out, file);
	strcat(file_out, add);
	if ((out = fopen(file_out, "wt")) == NULL)
	{
		printf("Ouput file can't be opened!\n");
		exit(1);
	}
	fprintf(out, "\t//%s\n", file);
	fprintf(out, "%d\t%f\t", size, fit);
	fprintf(out, "\n");
	for (i = 0; i < size; i++)
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
void town::fprint_allfi_mat(char* file, char* add, char* name, int len, double sga_min, double sga_raz, double* buf)
{
	int i;
	FILE* out;
	char file_out[500];
	strcpy(file_out, file);
	strcat(file_out, add);

	if ((out = fopen(file_out, "wt")) == NULL)
	{
		fprintf(out, "Ouput file can't be opened!\n");
		exit(1);
	}
	fprintf(out, "%s\n", name);
	fprintf(out, "%d\tLPD count\n", size);
	fprintf(out, "%d\tModel length\n", len);
	fprintf(out, "%.12f\tMinimum\n", sga_min);
	fprintf(out, "%.12f\tRazmah\n", sga_raz);
	for (i = 0; i < size; i++)
	{
		fprintf(out, "%d\t%d\t%.12f\t", tot[i].sta, tot[i].end, buf[i]);
		fprintf(out, "%d\t%s\n", tot[i].num, s[tot[i].num].oli);
	}
	fclose(out);
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
void DelHole(char* str)
{
	char* hole;
	hole = strstr(str, "\n");
	if (hole != NULL) *hole = 0;
	hole = strstr(str, "\r");
	if (hole != NULL) *hole = 0;
}
char* TransStr(char* d)
{
	int i, c, lens;
	lens = (int)strlen(d);
	for (i = 0; i < lens; i++)
	{
		c = int(d[i]);
		if (c < 97) d[i] = char(c + 32);
		//else break;
	}
	return(d);
}
int CheckStr(char* file, char* d, int n, int print, FILE* outlog)
{
	int i, len, ret;
	len = (int)strlen(d);
	ret = 1;
	for (i = 0; i < len; i++)
	{
		int di = (int)d[i];
		if (strchr("atgcATGC", di) != NULL)continue;
		if (strchr("nN", di) != NULL)
		{
			ret = 0; continue;
		}
		if (print == 1)fprintf(outlog, "File %s; sequence %d position %d (%c) bad. Sequence deleted!\n", file, n, i + 1, d[i]);
		ret = -1;
		break;
	}
	return(ret);
}
int IdeLet(char c)
{
	int ret;
	switch (c) {
	case 'a': ret = 0; break;
	case 'c': ret = 1; break;
	case 'g': ret = 2; break;
	case 't': ret = 3; break;
	case 'n': ret = -1; break;
	default: ret = -2;
	}
	return(ret);
}
void MinusStr1(int size, int j, int j0, double buf, double b[POPSIZE][POPSIZE])
{
	int k;
	b[j][j0] = 0;
	for (k = j0 + 1; k < size; k++)
	{
		if (b[j0][k] == 0)continue;
		b[j][k] -= buf * b[j0][k];
	}
}
void MinusStr2(int size, int j, int j0, double buf, double b[POPSIZE][POPSIZE])
{
	int k;
	b[j][j0] = 0;
	for (k = j - 1; k >= j0; k--)
	{
		if (b[j0][k] == 0)continue;
		b[j][k] -= buf * b[j0][k];
	}
}
void MinusStr(int size, int j, int j0, double buf, double b[POPSIZE][POPSIZE])
{
	int k;
	for (k = 0; k < size; k++)
	{
		if (b[j0][k] == 0)continue;
		b[j][k] -= buf * b[j0][k];
	}
}
int BackMat(int size, FILE* outlog)
{
	int i, j;
	double buf, b[POPSIZE][POPSIZE];
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
				fprintf(outlog, "\n1back %g(%d,%d) bij%g(%d,%d)\t", uw[i][j], i, j, b[i][j], i, j);
				//fprintf(outlog,"%g(%d)\n",uw[j][j],j);
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
				fprintf(outlog, "\n2back %g(%d,%d) bij%g(%d,%d)\t", uw[i][j], i, j, b[i][j], i, j);
				//fprintf(outlog,"%g(%d)\n",uw[j][j],j);
				//
				return -1;
			}
			uw[i][j] = 0;
		}
	}
	for (i = 0; i < size; i++)
	{
		buf = uw[i][i];
		if (fabs(buf) < 1e-030)
		{
			fprintf(outlog, "\nback3 %g %d\t", uw[i][i], i);
			return -1;
		}
		for (j = 0; j < size; j++)
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
	for (i = 0; i < size; i++)
	{
		for (j = 0; j < size; j++)
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
double EvalMahFIT(town* a, int n_train, int octa, int infc, int*** seq, int olen, double** dav, double** dcv, double** octa_prow, double* pexp, int* len, FILE* outlog)//double *hoxa_wei, qbs *qps, double **frp, double *octa_rat,
{
	int k, n, m;
	double f1[POPSIZE], df[POPSIZE], f2[POPSIZE], pfm[MOTLEN][16];

	int olen1 = olen - 1;
	for (k = 0; k < olen1; k++)for (n = 0; n < 16; n++)pfm[k][n] = 0;

	//printf("Enter Fit1 %d\n",n_train);
	//exit(1);
	for (k = 0; k < a->size; k++)
	{
		f2[k] = dav[a->tot[k].end - a->tot[k].sta][a->tot[k].num];
	}
	for (k = 0; k < a->size; k++)
	{
		for (n = 0; n < a->size; n++)uw[k][n] = 0;
	}
	for (k = 0; k < a->size; k++)f1[k] = 0;
	//printf("Enter Fit2 %d\n", n_train);
	for (m = 0; m < n_train; m++)
	{
		//if (len[m] < olen)continue;
		int ori = a->ori[m];
		int pos = a->pos[m];
		for (k = 0; k < olen1; k++)
		{
			int let = seq[ori][m][k + pos];
			if (let != -1)pfm[k][let]++;
		}
		double fs[POPSIZE];
		//printf("Seq %d\t",b);
		for (k = 0; k < a->size; k++)
		{
			int rlenk = (a->tot[k].end - a->tot[k].sta + 1);
			fs[k] = 0;
			for (n = a->tot[k].sta; n <= a->tot[k].end; n++)
			{
				int let = seq[ori][m][n + pos];
				if (let == -1)
				{
					a->fpr = 1;
					a->mah = a->fit = 0;
					return 0;
				}
				if (a->tot[k].num == let)fs[k]++;
			}
			fs[k] /= rlenk;
			//frp[b][k] = fs[k];
			//printf("%d %g\t", k, frp[b][k]);
		}
		//printf("\n");		
		for (k = 0; k < a->size; k++)
		{
			for (n = 0; n < a->size; n++)
			{
				uw[k][n] += fs[k] * fs[n];
			}
			f1[k] += fs[k];
		}
	}
	for (k = 0; k < olen1; k++)
	{
		for (n = 0; n < 16; n++)
		{
			pfm[k][n] /= n_train;
		}
	}
	if (infc == 0)a->inf = 1;
	else
	{
		a->inf = 0;
		for (k = 0; k < olen1; k++)
		{
			for (n = 0; n < 16; n++)
			{
				if (pfm[k][n] > 0)a->inf += pfm[k][n] * log(pfm[k][n] / pexp[n]);
			}
		}
		if (a->inf < 1)a->inf = 1;
		else a->inf = 1 + log10(a->inf);
	}
	for (k = 0; k < a->size; k++)for (n = 0; n < a->size; n++)uw[k][n] /= n_train;
	for (k = 0; k < a->size; k++)f1[k] /= n_train;
	for (k = 0; k < a->size; k++)
	{
		for (n = 0; n < a->size; n++)
		{
			uw[k][n] -= f1[k] * f1[n];
		}
	}
	for (k = 0; k < a->size; k++)
	{
		uw[k][k] += dcv[a->tot[k].end - a->tot[k].sta][a->tot[k].num];
	}
	for (k = 0; k < a->size; k++)for (n = 0; n < a->size; n++)uw[k][n] /= 2;
	//printf("Av1,2\n");
	for (k = 0; k < a->size; k++)
	{
		//av[k] = (f1[k] + f2[k]) / 2;
		df[k] = f1[k] - f2[k];
		//printf("%g %g\t", av[k], df[k]);
	}
	if (BackMat(a->size, outlog) == -1) { a->fit = 0; return 0; }
	a->mah = 0;
	for (k = 0; k < a->size; k++)
	{
		double buf = 0;
		for (n = 0; n < a->size; n++)buf += uw[k][n] * df[n];
		a->mah += buf * df[k];
	}
	//exit(1);
	//printf("Mah %f\n", a->mah);
	/*double c0 = 0;
	for (k = 0; k < a->size; k++)
	{
		buf[k] /= a->mah;
		c0 -= av[k] * buf[k];
		//printf("%d %g\t", k, buf[k]);
	}*/
	//printf("\nc0 %f\n", c0);
	//printf("Scores\n");//frp = NULL;
	//exit(1);
	/*for (b = 0; b < n_train; b++)//train
	{
		//printf("Seq %d\t", b);
		double sco = c0;
		for (k = 0; k < a->size; k++)
		{
			sco += buf[k] * frp[b][k];
		}
		qps[b].q = 1 - fabs(1 - sco);
		qps[b].n = xporti[b];
		//printf("%d %g\n", qps[b].n, qps[b].q);
	}*/
	//printf("\nScoresSort\n");
	//exit(1);
	//qsort(qps, n_train, sizeof(qps[0]), compare_qbs);
	/*for (b = 0; b < n_train; b++)
	{
		m = qps[b].n;
		printf("sco %f peak %d cep %d pos %d\n", qps[b].q,qps[b].n, a->ori[m], a->pos[m]);
	}*/
	double oct = 0;
	int sum = 0;
	//printf("ScoresWei\n");
	//exit(1);	
	//int olen1 = olen - 1;
	for (m = 0; m < n_train; m++)
	{
		//m = qps[b].n;
		//m = xporti[b];
		//	oct += octa_prow[a->ori[m]][m][a->pos[m]] * hoxa_wei[b];
		int posdi;
		if (a->ori[m] == 0)posdi = a->pos[m];
		else posdi = len[m] - olen - a->pos[m];
		double wei = octa_prow[m][posdi];
		if ((posdi < 0 || posdi > len[m] - olen) || (wei < -100 || wei > 100))// || posdi > len[m] - olen1)
		{
			fprintf(outlog, "ComplHoxaError peak %d ori %d len %d pos %d posdi %d hoxa %g\n", m, a->ori[m], a->pos[m], len[m], posdi, wei);
			exit(1);
		}
		/*if (wei == 0)
		{
			a->fpr = 1;
			a->mah = a->fit = 0;
			return 0;
		}*/
		oct += wei;
		sum++;
		//printf("oct %f sco %f peak %d cep %d pos %d pro %f wei %g\n", oct, qps[b].q, m, a->ori[m], a->pos[m], octa_prow[a->ori[m]][m][a->pos[m]], hoxa_wei[b]);
	}
	oct /= sum;
	//printf("Pow\n");
	oct = pow(10, oct);
	a->fpr = oct;
	a->fit = a->mah * a->fpr;
	//exit(1);
	return a->fit;
}
double EvalMahFITTrain(town* a, int nseq, int*** seq, char* file, int olen, int reg_max, int* len, double** dav, double** dcv, FILE* outlog, char* file_base, char* extmap, char* extmat, char *name)
{
	int k, n, m;
	double ret = 0;
	double av[POPSIZE], f1[POPSIZE], df[POPSIZE], f2[POPSIZE], buf[POPSIZE];

	for (k = 0; k < a->size; k++)
	{
		for (n = 0; n < a->size; n++)uw[k][n] = 0;
	}
	for (k = 0; k < a->size; k++)df[k] = 0;
	for (k = 0; k < a->size; k++)
	{
		f2[k] = dav[a->tot[k].end - a->tot[k].sta][a->tot[k].num];
	}
	for (k = 0; k < a->size; k++)f1[k] = 0;
	for (m = 0; m < nseq; m++)
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
		}
		for (k = 0; k < a->size; k++)
		{
			for (n = 0; n < a->size; n++)
			{
				uw[k][n] += fs[k] * fs[n];
			}
			f1[k] += fs[k];
		}
	}
	for (k = 0; k < a->size; k++)for (n = 0; n < a->size; n++)uw[k][n] /= nseq;
	for (k = 0; k < a->size; k++)df[k] /= nseq;
	for (k = 0; k < a->size; k++)
	{
		for (n = 0; n < a->size; n++)
		{
			uw[k][n] -= f1[k] * f1[n];
		}
	}
	for (k = 0; k < a->size; k++)
	{
		uw[k][k] += dcv[a->tot[k].end - a->tot[k].sta][a->tot[k].num];
	}
	for (k = 0; k < a->size; k++)for (n = 0; n < a->size; n++)uw[k][n] /= 2;
	{
		for (k = 0; k < a->size; k++)
		{
			av[k] = (f1[k] + f2[k]) / 2;
			df[k] = f1[k] - f2[k];
		}
	}
	if (BackMat(a->size, outlog) == -1) { a->fit = 0; return 0; }
	a->fit = 0;
	//double euc_buf=0;
	for (k = 0; k < a->size; k++)
	{
		buf[k] = 0;
		for (n = 0; n < a->size; n++)buf[k] += uw[k][n] * df[n];
		a->fit += buf[k] * df[k];
		//	euc_buf+=buf[k]*buf[k];
		//printf("%f\t%f\t%f\n",a->fit,buf,df[k]);
	}
	//	euc_buf=sqrt(euc_buf);
	//a->fit/=euc_buf;
	//double c0 = 0;
	for (k = 0; k < a->size; k++)
	{
		buf[k] /= a->fit;
		//c0 -= av[k] * buf[k];
	}	
	int b, o;
	double sga_min = 0, sga_max = 0;
	for (k = 0; k < a->size; k++)
	{
		if (buf[k] < 0)sga_min += buf[k];
		else sga_max += buf[k];
	}
	double sga_raz = sga_max - sga_min;
	a->fprint_allfi_mat(file_base, extmat, name, olen, sga_min, sga_raz, buf);
	char file_out[500];
	strcpy(file_out, file_base);
	strcat(file_out, extmap);
	FILE* out;
	if ((out = fopen(file_out, "wt")) == NULL)
	{
		printf("Ouput file can't be opened!\n");
		exit(1);
	}
	int olen1 = olen - 1;
	for (b = 0; b < nseq; b++)
	{
		double fs[POPSIZE];
		int lenp = len[b] - olen1;
		double sco_pos = 0;
		int obest = 0, pbest = 0;
		for (m = 0; m < lenp; m++)
		{			
			for (o = 0; o < 2; o++)
			{
				double sco = 0;
				for (k = 0; k < a->size; k++)
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
				sco = (sco - sga_min) / sga_raz;				
				if (sco > sco_pos)
				{
					sco_pos = sco;
					obest = o;
					pbest = m;
				}				
			}
		}		
		{
			char dir[] = "+-";
			char d[50];
			char letter[5] = "acgt";
			char err = 'n';
			memset(d, '\0', olen + 1);
			for (m = 0; m < olen1; m++)
			{
				int cod2 = seq[obest][b][pbest + m];
				int cod1 = cod2 / 4;
				if (cod1 >= 0 && cod1 <= 3)d[m] = letter[cod1];
				else d[m] = err;
			}
			{
				int cod2 = seq[obest][b][pbest + olen1 - 1];
				int cod1 = cod2 % 4;
				if (cod1 >= 0 && cod1 <= 3)d[olen1] = letter[cod1];
				else d[olen1] = err;
			}
			d[olen] = '\0';
			int pbest_sta = 0, pbest_end = 0;
			if (obest == 0)
			{
				pbest_sta = 1 + pbest;
				pbest_end = 1 + pbest + olen1;
			}
			else
			{
				pbest_sta = len[b] - pbest - olen1;
				pbest_end = len[b] - pbest;
			}

			fprintf(out, "%d\t%d\t%d\t%c\t%f\t%s\n", b + 1, pbest_sta, pbest_end, dir[obest], sco_pos, d);
		}
	}
	fclose(out);
	return a->fit;
}
int MutOlig0(town* a)
{
	int cy = 0, r1, num1, i, gom;
	do
	{
		gom = 0;
		r1 = rand() % a->size;
		num1 = rand() % 15;
		if (num1 >= a->tot[r1].num)num1++;
		int i0 = a->sum(num1);
		for (i = i0; i < i0 + a->deg[num1]; i++)
		{
			//if(a->tot[i].num!=num1)continue;
			//if(i==r1)continue;
			if (a->tot[r1].sta >= a->tot[i].sta && a->tot[r1].sta <= a->tot[i].end) { gom = 1; break; }
			if (a->tot[r1].end >= a->tot[i].sta && a->tot[r1].end <= a->tot[i].end) { gom = 1; break; }
			if (a->tot[i].sta >= a->tot[r1].sta && a->tot[i].sta <= a->tot[r1].end) { gom = 1; break; }
			if (a->tot[i].end >= a->tot[r1].sta && a->tot[i].end <= a->tot[r1].end) { gom = 1; break; }
		}
		cy++;
		if (cy > 10)return -1;
	} while (gom == 1);
	if (gom == 1)return -1;
	a->deg[a->tot[r1].num]--;
	a->deg[num1]++;
	a->tot[r1].num = num1;
	return r1;
}
/*
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
	for (j = 0; j < a->size; j++)
	{
		if (a->tot[j].sta > a->tot[j].end)
		{
			printf("LEN STA %d\tEND %d\n", a->tot[j].sta, a->tot[j].end);
			return -1;
		}
		if (a->tot[j].end - a->tot[j].sta + 1 > olen)
		{
			printf("LEN STA %d\tEND %d\n", a->tot[j].sta, a->tot[j].end);
			return -1;
		}
	}
	for (j = 0; j < a->size; j++)
	{

		if (a->tot[j].sta >= len || a->tot[j].sta < 0) { n[0] = j; return -1; }
		if (a->tot[j].end >= len || a->tot[j].end < 0) { n[0] = j; return -1; }
		odgs[a->tot[j].end - a->tot[j].sta]++;
	}
	for (i = 0; i < olen; i++)
	{
		if (a->odg[i] != odgs[i])
		{
			printf("Oli %d Announced %d Really %d", i + 1, a->odg[i], odgs[i]);
			return -1;
		}
	}
	for (j = 0; j < 16; j++)
	{
		if (a->deg[j] == 0)continue;
		//if(a->deg[j]==1){ns++;continue;}
		gom = 0;
		for (i = ns; i < ns + a->deg[j]; i++)
		{
			for (k = ns; k < ns + a->deg[j]; k++)
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
	for (i = 0; i < olen; i++)
	{
		if (a->odg[i] < 0)
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
}*/
int MutCry0(town* a, int oln, int reg_max)
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
			if (a->tot[n1m1].end < oln2)a1 = a->tot[n1m1].end + 1;
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
		if (cy > 10)return -1;
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
int MutRegShift(town* a, int nseq, int* len, int olen, int& npeak, int& nori, int& npos)
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
}

/*int MutRegShiftHoxa(town *a, int n_train, int *xporti, int &npeak, int &nori, int &npos, int *len_octa, int **octa_prowb, int *len, int olen)
{
	//printf("In Peak %d Ori %d Pos %d ",npeak,nori,npos);
	int r1, r2, oln1 = olen - 1;
	r1 = rand() % n_train;
	npeak = xporti[r1];
	int cep = a->ori[npeak], cep1 = 1 - cep;
	int inx = len_octa[npeak] + len_octa[npeak] - 1;
	r2 = rand() % inx;
	if (r2 < len_octa[npeak])
	{
		nori = cep1;
		npos = octa_prowb[npeak][r2];
	}
	else
	{
		nori = cep;
		r2 -= len_octa[npeak];
		int r3 = octa_prowb[npeak][r2];
		if (r3 >= a->pos[npeak])npos = octa_prowb[npeak][r2 + 1];
		else npos = r3;
	}
	a->ori[npeak] = nori;
	a->pos[npeak] = npos;
	//a->posc[npeak] = len[npeak] - olen - npos;
	//printf("Out2 Peak %d Pos %d Ori %d",npeak, npos,nori);
	return 1;

}*/
int MutRegShiftHoxaW(town* a, int nseq, int& npeak, int& nori, int& npos, int* len_octa, int** octa_prowb, int** octa_prows)
{
	//printf("In Peak %d Ori %d Pos %d ",npeak,nori,npos);
	int r2, j;
	npeak = rand() % nseq;
	int inx = octa_prows[npeak][len_octa[npeak] - 1];
	r2 = rand() % inx;
	for (j = 0; j < len_octa[npeak]; j++)
	{
		if (octa_prows[npeak][j] > r2)
		{
			npos = octa_prowb[npeak][j];
			break;
		}
	}
	if (npos == a->pos[npeak])nori = 1 - a->ori[npeak];
	else nori = rand() % 2;
	a->ori[npeak] = nori;
	a->pos[npeak] = npos;
	//a->posc[npeak] = len[npeak] - olen - npos;
	//printf("Out2 Peak %d Pos %d Ori %d",npeak, npos,nori);
	return 1;
}
int RecFeat(town a1, town a2, int(*cop)[2], int max)
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
			for (i = 0; i < ret; i++)
			{
				if (r1 >= cop[i][0])r1++;
			}
			num1 = a1.tot[r1].num;
			if (a2.deg[num1] == 0)break;
			int i0 = a2.sum(num1);
			for (i = i0; i < i0 + a2.deg[num1]; i++)
			{
				if (gom1 > 1)break;
				if (a2.tot[i].sta > a1.tot[r1].end)break;
				if (a2.tot[i].end < a1.tot[r1].sta)continue;
				if (a2.tot[i].sta >= a1.tot[r1].sta && a2.tot[i].sta <= a1.tot[r1].end) { gom1++; r0 = i; continue; }
				if (a2.tot[i].end >= a1.tot[r1].sta && a2.tot[i].end <= a1.tot[r1].end) { gom1++; r0 = i; continue; }
				if (a1.tot[r1].sta >= a2.tot[i].sta && a1.tot[r1].sta <= a2.tot[i].end) { gom1++; r0 = i; continue; }
				if (a1.tot[r1].end >= a2.tot[i].sta && a1.tot[r1].end <= a2.tot[i].end) { gom1++; r0 = i; continue; }
			}
			if (gom1 == 1)
			{
				if (a1.tot[r1].end == a2.tot[r0].end && a1.tot[r1].sta == a2.tot[r0].sta)
				{
					gom1++;
				}
			}
			if (cy > reg[0])return ret;
			cy++;
		} while (gom1 > 1);
		cy = 0;
		do
		{
			gom2 = 0;
			if (gom1 == 0)r2 = rand() % (a2.size - ret);
			else r2 = r0;
			for (i = 0; i < ret; i++)
			{
				if (r2 >= cop[i][1])r2++;
			}
			num2 = a2.tot[r2].num;
			if (a1.deg[num2] == 0)break;
			if (num1 == num2 && (a1.tot[r1].end == a2.tot[r2].end && a1.tot[r1].sta == a2.tot[r2].sta))continue;
			int i0 = a1.sum(num2);
			for (i = i0; i < i0 + a1.deg[num2]; i++)
			{
				if (gom2 > 1)break;
				if (a1.tot[i].sta > a2.tot[r2].end)break;
				if (a1.tot[i].end < a2.tot[r2].sta)continue;
				if (a1.tot[i].sta >= a2.tot[r2].sta && a1.tot[i].sta <= a2.tot[r2].end) { gom2++; continue; }
				if (a1.tot[i].end >= a2.tot[r2].sta && a1.tot[i].end <= a2.tot[r2].end) { gom2++; continue; }
				if (a2.tot[r2].sta >= a1.tot[i].sta && a2.tot[r2].sta <= a1.tot[i].end) { gom2++; continue; }
				if (a2.tot[r2].end >= a1.tot[i].sta && a2.tot[r2].end <= a1.tot[i].end) { gom2++; continue; }
			}
			if (gom1 == 1 && gom2 == 1)break;
			if (cy > reg[1])return ret;
			cy++;
		} while (gom2 != 0);
		cop[ret][0] = r1;
		cop[ret][1] = r2;
		//		a1.tot[r1].print_all();
		//		a2.tot[r2].print_all();
		//		printf("\n");
		ret++;
	} while (ret < max);
	//	a1.print_all();
	//	a2.print_all();
	return ret;
}
//original: simple exchange
int Reco2_Original(town* a1, town* a2, int* nsi, int* num)
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
		for (i = i0; i < i0 + a2->deg[num1]; i++)
		{
			if (gom1 > 1)break;
			if (a2->tot[i].sta > a1->tot[r1].end)break;
			if (a2->tot[i].end < a1->tot[r1].sta)continue;
			if (a2->tot[i].sta >= a1->tot[r1].sta && a2->tot[i].sta <= a1->tot[r1].end) { gom1++; r0 = i; continue; }
			if (a2->tot[i].end >= a1->tot[r1].sta && a2->tot[i].end <= a1->tot[r1].end) { gom1++; r0 = i; continue; }
			if (a1->tot[r1].sta >= a2->tot[i].sta && a1->tot[r1].sta <= a2->tot[i].end) { gom1++; r0 = i; continue; }
			if (a1->tot[r1].end >= a2->tot[i].sta && a1->tot[r1].end <= a2->tot[i].end) { gom1++; r0 = i; continue; }
		}
		if (gom1 == 1)
		{
			if (a1->tot[r1].end == a2->tot[r0].end && a1->tot[r1].sta == a2->tot[r0].sta)
			{
				gom1++;
			}
		}
		if (cy > 10)return -1;
		cy++;
	} while (gom1 > 1);
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
		for (i = i0; i < i0 + a1->deg[num2]; i++)
		{
			if (gom2 > 1)break;
			if (a1->tot[i].sta > a2->tot[r2].end)break;
			if (a1->tot[i].end < a2->tot[r2].sta)continue;
			if (a1->tot[i].sta >= a2->tot[r2].sta && a1->tot[i].sta <= a2->tot[r2].end) { gom2++; continue; }
			if (a1->tot[i].end >= a2->tot[r2].sta && a1->tot[i].end <= a2->tot[r2].end) { gom2++; continue; }
			if (a2->tot[r2].sta >= a1->tot[i].sta && a2->tot[r2].sta <= a1->tot[i].end) { gom2++; continue; }
			if (a2->tot[r2].end >= a1->tot[i].sta && a2->tot[r2].end <= a1->tot[i].end) { gom2++; continue; }
		}
		if (gom1 == 1 && gom2 == 1)break;
		if (cy > 10)return -1;
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
int Reco2_Economic(town* a1, town* a2, int* nsi, int* num)
{
	int i0, j, r1, r2, num1, num2, i, gom1 = 0, gom2 = 0, ret;
	int size1 = Min(a1->size, a2->size);
	int ord[POPSIZE];
	for (i = 0; i < size1; i++)ord[i] = i;
	BigMixI(ord, size1);
	size1 /= 2;
	for (j = 0; j < size1; j++)
	{
		r1 = ord[j];
		gom1 = 0;
		num1 = a1->tot[r1].num;
		if (a2->deg[num1] == 0)continue;
		i0 = a2->sum(num1);
		if (a2->deg[num1] == 1)r2 = i0;
		else
		{
			for (i = i0; i < i0 + a2->deg[num1]; i++)
			{
				if (a2->tot[i].end < a1->tot[r1].sta)continue;
				if (gom1 > 1)break;
				if (a2->tot[i].sta > a1->tot[r1].end)break;
				if (a2->tot[i].sta >= a1->tot[r1].sta && a2->tot[i].sta <= a1->tot[r1].end) { gom1++; r2 = i; continue; }//1s 2s 1e
				if (a2->tot[i].end >= a1->tot[r1].sta && a2->tot[i].end <= a1->tot[r1].end) { gom1++; r2 = i; continue; }//1s 2e 1e
				if (a1->tot[r1].sta >= a2->tot[i].sta && a1->tot[r1].sta <= a2->tot[i].end) { gom1++; r2 = i; continue; }//2s 1s 2e
				if (a1->tot[r1].end >= a2->tot[i].sta && a1->tot[r1].end <= a2->tot[i].end) { gom1++; r2 = i; continue; }//2s 1e 2e
			}
		}
		if (gom1 > 1)continue;
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
		for (i = i0; i < i0 + a1->deg[num2]; i++)
		{
			if (i == r1)continue;
			if (gom2 > 1)break;
			if (a1->tot[i].sta > a2->tot[r2].end)break;
			if (a1->tot[i].end < a2->tot[r2].sta)continue;
			if (a1->tot[i].sta >= a2->tot[r2].sta && a1->tot[i].sta <= a2->tot[r2].end) { gom2++; continue; }
			if (a1->tot[i].end >= a2->tot[r2].sta && a1->tot[i].end <= a2->tot[r2].end) { gom2++; continue; }
			if (a2->tot[r2].sta >= a1->tot[i].sta && a2->tot[r2].sta <= a1->tot[i].end) { gom2++; continue; }
			if (a2->tot[r2].end >= a1->tot[i].sta && a2->tot[r2].end <= a1->tot[i].end) { gom2++; continue; }
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
int Reco2_One_dinucleotide_full(town* a1, town* a2)
{
	int gom = -1, i, j, w1, w2, i1, i2;
	int size1 = 16;
	int ord[16];
	for (i = 0; i < size1; i++)ord[i] = i;
	BigMixI(ord, size1);
	for (j = 0; j < size1; j++)
	{
		if (a1->deg[ord[j]] == 0)continue;
		if (a1->deg[ord[j]] != a2->deg[ord[j]])continue;
		i1 = a1->sum(ord[j]);
		i2 = a2->sum(ord[j]);
		for (i = 0; i < a1->deg[ord[j]]; i++)
		{
			int k1 = i1 + i;
			int k2 = i2 + i;
			if (a1->tot[k1].sta != a2->tot[k2].sta) { gom = ord[j]; break; }
			if (a1->tot[k1].end != a2->tot[k2].end) { gom = ord[j]; break; }
		}
		if (gom != -1)break;
	}
	if (gom == -1)return -1;
	uno buf;
	for (j = 0; j < a1->deg[gom]; j++)
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
//total: exchange of subsets of two dinucl.types
int Reco2_Two_dinucleotides_full(town* a1, town* a2, int nseq, int reg_max)
{
	int i, j, k, gom = 0, x = -1, y = -1;
	int ord[16];
	for (i = 0; i < 16; i++)ord[i] = i;
	BigMixI(ord, 16);
	for (j = 0; j < 15; j++)
	{
		x = ord[j];
		if (a1->deg[x] == 0 || a2->deg[x] == 0)continue;
		int goma = 0;
		if (a1->deg[x] == a2->deg[x])
		{
			int sj1 = a1->sum(x), sj2 = a2->sum(x);
			for (k = 0; k < a1->deg[x]; k++)
			{
				int k1 = sj1 + k;
				int k2 = sj2 + k;
				if (a1->tot[k1].sta != a2->tot[k2].sta) { goma = 1; break; }
				if (a1->tot[k1].end != a2->tot[k2].end) { goma = 1; break; }
			}
		}
		else goma = 1;
		if (goma == 0)continue;
		int d1 = a1->deg[x] - a2->deg[x];
		for (i = j + 1; i < 16; i++)
		{
			y = ord[i];
			if (a1->deg[y] == 0 || a2->deg[y] == 0)continue;
			int d2 = a2->deg[y] - a1->deg[y];
			if (d2 != d1)continue;//if (a1->deg[x] + a1->deg[y] != a2->deg[x] + a2->deg[y])continue;
			int gomb = 0;
			if (a1->deg[y] == a2->deg[y])
			{
				int si1 = a1->sum(y), si2 = a2->sum(y);
				for (k = 0; k < a1->deg[y]; k++)
				{
					int k1 = si1 + k;
					int k2 = si2 + k;
					if (a1->tot[k1].sta != a2->tot[k2].sta) { gomb = 1; break; }
					if (a1->tot[k1].end != a2->tot[k2].end) { gomb = 1; break; }
				}
			}
			else gomb = 1;
			if (gomb == 1)
			{
				gom = 1;
				break;
			}
		}
		if (gom == 1)
		{
			break;
		}
	}
	if (gom == 0)return -1;
	else
	{
		town b1, b2;	//old ones	
		for (i = 0; i < 16; i++)
		{
			b1.deg[i] = a1->deg[i];
			b2.deg[i] = a2->deg[i];
		}
		for (i = 0; i < a1->size; i++)
		{
			b1.tot[i].sta = a1->tot[i].sta;
			b1.tot[i].end = a1->tot[i].end;
			b1.tot[i].num = a1->tot[i].num;
			//a1->tot[i].get_copy(&b1.tot[i]);
		}
		for (i = 0; i < a2->size; i++)
		{
			b2.tot[i].sta = a2->tot[i].sta;
			b2.tot[i].end = a2->tot[i].end;
			b2.tot[i].num = a2->tot[i].num;
			//	a2->tot[i].get_copy(&b2.tot[i]);
		}
		for (j = 0; j < a1->deg[x]; j++)
		{
			int w = a1->tot[j].end - a1->tot[j].sta;
			a1->odg[w]--;
		}
		for (j = 0; j < a2->deg[x]; j++)
		{
			int w = a2->tot[j].end - a2->tot[j].sta;
			a2->odg[w]--;
		}
		for (j = 0; j < a1->deg[y]; j++)
		{
			int w = a1->tot[j].end - a1->tot[j].sta;
			a1->odg[w]--;
		}
		for (j = 0; j < a2->deg[y]; j++)
		{
			int w = a2->tot[j].end - a2->tot[j].sta;
			a2->odg[w]--;
		}
		int ka1 = 0, ka2 = 0, kb1 = 0, kb2 = 0;
		for (i = 0; i < 16; i++)
		{
			if (i == x || i == y)
			{
				a1->deg[i] = b2.deg[i];
				a2->deg[i] = b1.deg[i];
				for (j = 0; j < a1->deg[i]; j++)
				{
					//b2.tot[kb2].get_copy(&a1->tot[ka1]);
					a1->tot[ka1].sta = b2.tot[kb2].sta;
					a1->tot[ka1].end = b2.tot[kb2].end;
					a1->tot[ka1].num = b2.tot[kb2].num;
					kb2++;
					ka1++;
				}
				for (j = 0; j < a2->deg[i]; j++)
				{
					//b1.tot[kb1].get_copy(&a2->tot[ka2]); 
					a2->tot[ka2].sta = b1.tot[kb1].sta;
					a2->tot[ka2].end = b1.tot[kb1].end;
					a2->tot[ka2].num = b1.tot[kb1].num;
					kb1++;
					ka2++;
				}
			}
			else
			{
				a1->deg[i] = b1.deg[i];
				a2->deg[i] = b2.deg[i];
				for (j = 0; j < a1->deg[i]; j++)
				{
					//b1.tot[kb1].get_copy(&a1->tot[ka1]); 
					a1->tot[ka1].sta = b1.tot[kb1].sta;
					a1->tot[ka1].end = b1.tot[kb1].end;
					a1->tot[ka1].num = b1.tot[kb1].num;
					kb1++;
					ka1++;
				}
				for (j = 0; j < a2->deg[i]; j++)
				{
					//b2.tot[kb2].get_copy(&a2->tot[ka2]); 
					a2->tot[ka2].sta = b2.tot[kb2].sta;
					a2->tot[ka2].end = b2.tot[kb2].end;
					a2->tot[ka2].num = b2.tot[kb2].num;
					kb2++;
					ka2++;
				}
			}
		}
		for (j = 0; j < a1->deg[x]; j++)
		{
			int w = a1->tot[j].end - a1->tot[j].sta;
			a1->odg[w]++;
		}
		for (j = 0; j < a2->deg[x]; j++)
		{
			int w = a2->tot[j].end - a2->tot[j].sta;
			a2->odg[w]++;
		}
		for (j = 0; j < a1->deg[y]; j++)
		{
			int w = a1->tot[j].end - a1->tot[j].sta;
			a1->odg[w]++;
		}
		for (j = 0; j < a2->deg[y]; j++)
		{
			int w = a2->tot[j].end - a2->tot[j].sta;
			a2->odg[w]++;
		}
		return 1;
	}
}
//canonical: : classic recombination within one dinucl.type
int Reco2_One_dinucleotide_local(town* a1, town* a2, int olen, int reg_max)
{
	int i, j, j0, w1, w2, i1, i2;
	int size1 = 16;
	int ord[16];
	int rpt = 0, rp[SEQLEN];
	for (i = 0; i < size1; i++)ord[i] = i;
	BigMixI(ord, size1);
	for (j = 0; j < size1; j++)
	{
		if (a1->deg[ord[j]] <= 1)continue;
		if (a1->deg[ord[j]] != a2->deg[ord[j]])continue;
		i1 = a1->sum(ord[j]);
		i2 = a2->sum(ord[j]);
		for (i = 0; i < a1->deg[ord[j]]; i++)
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
					if (over < 0)continue;
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
		if (w1 < 0 || w1 >= reg_max)return -1;
		w2 = a1->tot[k1].end - a2->tot[k2].sta;
		if (w2 < 0 || w2 >= reg_max)return -1;
		a1->odg[w1]++;
		a2->odg[w2]++;
		int end1 = a1->tot[k1].end;
		a1->tot[k1].end = a2->tot[k2].end;
		a2->tot[k2].end = end1;
	}
	uno buf;
	int cha = 0;
	for (j = rp[rpt0] + 1; j < a1->deg[j0]; j++)
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
int Reco2Peak(town* a1, town* a2, int n_train)
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

void MixPop(town* a, town* b)
{
	town c = *a;
	*a = *b;
	*b = c;
}
int GomTown(town a, town b, int nseq, int check_lpd)
{
	int i;
	for (i = 0; i < nseq; i++)
	{
		if (b.pos[i] != a.pos[i])return 0;
		if (b.ori[i] != a.ori[i])return 0;
	}
	if (check_lpd == 0)return -1;
	for (i = 0; i < a.size; i++)
	{
		if (b.tot[i].sta != a.tot[i].sta)return 0;
		if (b.tot[i].end != a.tot[i].end)return 0;
		if (b.tot[i].num != a.tot[i].num)return 0;
	}
	return -1;
}
int GomTown2(town a, town b)
{
	int i;
	for (i = 0; i < a.size; i++)
	{
		if (b.tot[i].sta != a.tot[i].sta)return 0;
		if (b.tot[i].end != a.tot[i].end)return 0;
		if (b.tot[i].num != a.tot[i].num)return 0;
	}
	return -1;
}
void DelChar(char* str, char c)
{
	int i, lens, size;

	size = 0;
	lens = (int)strlen(str);
	for (i = 0; i < lens; i++)
	{
		if (str[i] != c)str[size++] = str[i];
	}
	str[size] = '\0';
}
int ComplStr(char* d)
{
	int i, len;
	len = (int)strlen(d);
	char d1[SEQLEN];
	strcpy(d1, d);
	//	memset(d,0,sizeof(d));
	for (i = 0; i < len; i++)
	{
		switch (d1[len - i - 1])
		{
		case 'a': { d[i] = 't'; break; }
		case 't': { d[i] = 'a'; break; }
		case 'c': { d[i] = 'g'; break; }
		case 'g': { d[i] = 'c'; break; }
		default: d[i] = 'n';
		}
	}
	return 1;
}
void GetSost(char* d, int word, int size, int* sost)
{
	int i, j, k, i_sost, let;
	char letter[5] = "acgt";
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
void GetSostPro(char* d, int word, int size, int* sost, double* sost_pro)
{
	int i, j, k, i_sost;
	char letter[6] = "acgtn";
	int ten[8] = { 1, 4, 16, 64, 256, 1024, 4096, 16384 };
	int lens = (int)strlen(d);
	int word1 = word - 1;
	for (i = 0; i < size; i++)sost[i] = 0;
	for (i = 0; i < lens - word1; i++)
	{
		i_sost = 0;
		int err = 1;
		for (j = word1; j >= 0; j--)
		{
			for (k = 0; k < 4; k++)
			{
				if (d[i + j] == letter[k])
				{
					i_sost += ten[word1 - j] * k;
					err = 0;
					break;
				}
			}
			if (err == 1)break;
		}
		if (err == 0)
		{
			sost[i_sost]++;
			sost_pro[i] = (double)i_sost;
		}
		else sost_pro[i] = -1;
	}
}
void EvalSeq(char* file, int& nseq, int olen, int len_peak_max, FILE* outlog)
{
	char l[SEQLEN], d[SEQLEN], head[400];
	int fl = 0;
	FILE* in;

	if ((in = fopen(file, "rt")) == NULL)
	{
		printf("EvalSeq! Input file %s can't be opened!\n", file);
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
			int lenx = (int)strlen(d);
			int check = CheckStr(file, d, n, 1, outlog);
			if ((lenx >= olen && lenx <= len_peak_max) && check != -1)nseq++;
			if (fl == -1)
			{
				fclose(in);
				break;
			}
			n++;
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
void EvalLen(char* file, int* len, int olen, int len_peak_max, FILE* outlog)
{
	char l[SEQLEN], d[SEQLEN], head[400];
	int fl = 0;
	FILE* in;

	if ((in = fopen(file, "rt")) == NULL)
	{
		printf("EvalLen! Input file %s can't be opened!\n", file);
		exit(1);
	}
	char symbol = fgetc(in);
	rewind(in);
	int nn = 0;
	int n = 0;
	while (n >= 0)
	{
		if (fgets(l, sizeof(l), in) == NULL) fl = -1;
		if (*l == '\n' && fl != -1)continue;
		if (((*l == symbol) || (fl == -1)) && (fl != 0))
		{
			int lenx = (int)strlen(d);
			int check = CheckStr(file, d, nn, 0, outlog);
			if ((lenx >= olen && lenx <= len_peak_max) && check != -1)len[n++] = lenx;
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
void ReadSeq(char* file, int nseq, int* len, int*** seq_real, int olen, double* octaf, int* octa1, int octa, int octa_size, double** octa_pro1, int len_peak_max, FILE* outlog)
{
	char l[SEQLEN], d[2][SEQLEN], head[400];
	int fl = 0, i, j;
	FILE* in;

	if ((in = fopen(file, "rt")) == NULL)
	{
		printf("ReadSeq! Input file %s can't be opened!\n", file);
		exit(1);
	}
	char symbol = fgetc(in);
	rewind(in);
	int size = 1;
	for (i = 0; i < octa; i++)size *= 4;

	int nn = 0, n = 0;
	while (n >= 0)
	{
		if (fgets(l, sizeof(l), in) == NULL) fl = -1;
		if (*l == '\n' && fl != -1)continue;
		if (((*l == symbol) || (fl == -1)) && (fl != 0))
		{
			int lenx = (int)strlen(d[0]);
			int check = CheckStr(file, d[0], nn, 0, outlog);
			nn++;
			if ((lenx >= olen && lenx <= len_peak_max) && check != -1)
			{
				TransStr(d[0]);
				d[0][len[n]] = '\0';
				strcpy(d[1], d[0]);
				ComplStr(d[1]);
				d[1][len[n]] = '\0';
				GetSostPro(d[0], octa, octa_size, octa1, octa_pro1[n]);
				for (i = 0; i < octa_size; i++)octaf[i] += octa1[i];
				GetSost(d[1], octa, octa_size, octa1);
				for (i = 0; i < octa_size; i++)octaf[i] += octa1[i];
				for (j = 0; j < 2; j++)
				{
					int cod[2];
					cod[0] = IdeLet(d[j][0]);
					int len1 = len[n] - 1;
					for (i = 0; i < len1; i++)
					{
						cod[1] = IdeLet(d[j][i + 1]);
						if (cod[0] >= 0 && cod[1] >= 0)seq_real[j][n][i] = 4 * cod[0] + cod[1];
						else
						{
							seq_real[j][n][i] = -1;
							//	printf("ReadSeq! Input file %s peak %d strand %d position %d error!\n%s\n", file, n, j, i, d[j]);
								//exit(1);
						}
						cod[0] = cod[1];
					}
				}
				n++;
			}
			else
			{
				if (lenx < olen)
				{
					printf("ReadSeq! Short peak %d (Len %d) ignored\n", nn + 1, lenx);
				}
				if (lenx > len_peak_max)
				{
					printf("ReadSeq! Long peak %d (Len %d) ignored\n", nn + 1, lenx);
				}
				if (check != 1)
				{
					printf("ReadSeq! Unusual symbol, peak %d partially ignored\n%s\n", nn + 1, d[0]);
				}
			}
			if (fl == -1)
			{
				fclose(in);
				double sumf = 0;
				for (i = 0; i < octa_size; i++)
				{
					if (octaf[i] == 0)octaf[i] = 1;
					sumf += octaf[i];
				}
				for (i = 0; i < octa_size; i++)octaf[i] /= sumf;
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
void ReadSeqBack(char* file, int nseq, int* len, int olen, int reg_max, double** dav, double** dcv, double* octaf, double * pexp, int* octa1, int octa, int octa_size, int len_max, int len_peak_max, FILE* outlog)
{
	char l[SEQLEN], d[SEQLEN], head[400];
	int m, i, j, k, v, fl = 0;
	int n_backr[LPDLEN];
	double freq[16];
	FILE* in;

	int* seq_back;
	seq_back = new int[len_max - 1];

	if ((in = fopen(file, "rt")) == NULL)
	{
		printf("ReadSeq! Input file %s can't be opened!\n", file);
		exit(1);
	}
	char symbol = fgetc(in);
	rewind(in);
	for (j = 0; j < reg_max; j++)n_backr[j] = 0;
	int nn = 0, n = 0;
	int all_pos = 0;
	while (n >= 0)
	{
		if (fgets(l, sizeof(l), in) == NULL) fl = -1;
		if (*l == '\n' && fl != -1)continue;
		if (((*l == symbol) || (fl == -1)) && (fl != 0))
		{
			int lenx = (int)strlen(d);
			int check = CheckStr(file, d, nn, 0, outlog);
			nn++;
			if ((lenx >= olen && lenx <= len_peak_max) && check != -1)
			{
				TransStr(d);
				d[len[n]] = '\0';
				int cod[2];
				for (m = 0; m < 2; m++)
				{
					if (m == 1)ComplStr(d);
					GetSost(d, octa, octa_size, octa1);
					for (i = 0; i < octa_size; i++)octaf[i] += octa1[i];
					cod[0] = IdeLet(d[0]);
					for (k = 0; k < len[n] - 1; k++)
					{
						cod[1] = IdeLet(d[k + 1]);
						int code = -1;
						if (cod[0] >= 0 && cod[1] >= 0) code = 4 * cod[0] + cod[1];
						//	else
						//	{
								//printf("ReadSeqBack! Input file %s peak %d strand %d position %d error!\n%s\n", file, n, m, k, d);
								//exit(1);
						//	}
						seq_back[k] = code;
						cod[0] = cod[1];
					}
					for (j = 1; j <= reg_max; j++)
					{
						int j1 = j - 1;
						for (k = 0; k < len[n] - j; k++)
						{
							int gom = 1;
							for (i = 0; i < 16; i++)freq[i] = 0;
							for (v = 0; v < j; v++)
							{
								int sb = seq_back[k + v];
								if (sb == -1)
								{
									gom = -1;
									break;
								}
								else
								{
									freq[sb]++;
									pexp[sb]++;
									all_pos++;
								}
							}
							if (gom == 1)
							{
								n_backr[j1]++;
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
				}
				n++;
			}
			else
			{
				if (lenx < olen)
				{
					printf("ReadSeqBack! Short peak %d (Len %d) ignored\n", nn + 1, lenx);
				}
				if (lenx > len_peak_max)
				{
					printf("ReadSeqBack! Long peak %d (Len %d) ignored\n", nn + 1, lenx);
				}
				if (check != 1)
				{
					printf("ReadSeqBack!Unusual symbol, peak %d partially ignored\n%s\n", nn + 1, d);
					//exit(1);
				}
			}
			if (fl == -1)
			{

				for (j = 0; j < reg_max; j++)
				{
					n_backr[j] *= 2;
					for (k = 0; k < 16; k++)
					{
						dav[j][k] /= n_backr[j];
						dcv[j][k] /= n_backr[j];
						dcv[j][k] -= dav[j][k] * dav[j][k];
					}
				}
				fclose(in);
				double sumf = 0;
				for (i = 0; i < octa_size; i++)
				{
					if (octaf[i] == 0)octaf[i] = 1;
					sumf += octaf[i];
				}
				for (i = 0; i < octa_size; i++)octaf[i] /= sumf;
				for (k = 0; k < 16; k++)pexp[k] /= all_pos;
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
int main(int argc, char* argv[])
{
	int* len, nseq, * lenb, nseqb, i, j, k, m, n;
	char file_for[500], file_back[500], path_fasta[500], path_out[500], pfile_for[500], pfile_back[500], file_log[500];
	char filef[500], file_out[500];
	int*** seq_real;	
	double** dav;//dinucl.content background
	double** dcv;//self covariations for regions LPD	
	//	double **frp;//LPD frequencies
	double* qp;//train scores	
	int** octa_prowb, * len_octa, ** octa_prows;// octa position lists, octa position counts
	double** octa_pro1, ** octa_prow, * thr_octa;
	double pexp[16];

	FILE* outlog;

	if (argc != 12)
	{
		puts("Sintax: 1path_both_fasta 2file_forground 3file_background 4int max_LPD_length 5int motif_len 6int size 7int olig_background 8int infc(0no,1yes) 9path_out 10int max_peak_len 11file output log");//  5<pop_size>
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
	int octa = atoi(argv[7]);
	int infc = atoi(argv[8]); // information content in GA fitness
	strcpy(path_out, argv[9]);
	int len_peak_max = atoi(argv[10]); //2500;
	strcpy(file_log, path_out);
	strcat(file_log, argv[11]);
	if ((outlog = fopen(file_log, "wt")) == NULL)
	{
		fprintf(outlog, "Input file %s can't be opened!\n", file_log);
		exit(1);
	}

	int olen1 = olen - 1;
	srand((unsigned)time(NULL));
	if (size > POPSIZE)
	{
		printf("Maximal population size too large %d\n", size);
		exit(1);
	}
	time_t tnow = time(NULL);
	fprintf(outlog, "%s", ctime(&tnow));
	dcv = new double* [reg_max];
	if (dcv == NULL) return -1;
	for (i = 0; i < reg_max; i++)
	{
		dcv[i] = new double[16];
		if (dcv[i] == NULL) return -1;
	}
	dav = new double* [reg_max];
	if (dav == NULL) return -1;
	for (i = 0; i < reg_max; i++)
	{
		dav[i] = new double[16];
		if (dav[i] == NULL) return -1;
	}
	//foreground
	EvalSeq(pfile_for, nseq, olen, len_peak_max, outlog);
	len = new int[nseq];
	if (len == NULL) { puts("Out of memory..."); exit(1); }
	EvalLen(pfile_for, len, olen, len_peak_max, outlog);
	//background
	EvalSeq(pfile_back, nseqb, olen, len_peak_max, outlog);
	lenb = new int[nseqb];
	if (lenb == NULL) { puts("Out of memory..."); exit(1); }
	EvalLen(pfile_back, lenb, olen, len_peak_max, outlog);
	seq_real = new int** [2];
	if (seq_real == NULL) { puts("Out of memory..."); exit(1); }
	for (i = 0; i < 2; i++)
	{
		seq_real[i] = new int* [nseq];
		for (j = 0; j < nseq; j++)
		{
			seq_real[i][j] = new int[len[j] - 1];
			if (seq_real[i][j] == NULL) { puts("Out of memory..."); exit(1); }
		}
	}
	for (i = 0; i < 2; i++)
	{
		for (j = 0; j < nseq; j++)for (k = 0; k < len[j] - 1; k++)seq_real[i][j][k] = -1;
	}
	for (j = 0; j < reg_max; j++)
	{
		for (i = 0; i < 16; i++)dcv[j][i] = dav[j][i] = 0;
	}
	int octa_size = 1;
	for (i = 0; i < octa; i++)octa_size *= 4;
	octa_pro1 = new double* [nseq];
	if (octa_pro1 == NULL) return -1;
	for (i = 0; i < nseq; i++)
	{
		octa_pro1[i] = new double[len[i] - octa + 1];
		if (octa_pro1[i] == NULL) return -1;
	}
	octa_prow = new double* [nseq];
	if (octa_prow == NULL) return -1;
	for (i = 0; i < nseq; i++)
	{
		octa_prow[i] = new double[len[i] - octa + 1];
		if (octa_prow[i] == NULL) return -1;
	}
	int len_max = 0;
	for (i = 0; i < nseq; i++)if (len[i] > len_max)len_max = len[i];
	double* octa_pro1p;
	octa_pro1p = new double[len_max];
	if (octa_pro1p == NULL) { puts("Out of memory..."); exit(1); }
	fprintf(outlog, "Start1\n");
	double* octa_rat;
	octa_rat = new double[octa_size];
	if (octa_rat == NULL) { puts("Out of memory..."); exit(1); }
	double* octa_av;
	octa_av = new double[octa_size];
	if (octa_av == NULL) { puts("Out of memory..."); exit(1); }
	int* octa1;
	octa1 = new int[octa_size];
	if (octa1 == NULL) { puts("Out of memory..."); exit(1); }
	for (i = 0; i < octa_size; i++)octa_av[i] = 0;
	for (i = 0; i < octa_size; i++)octa1[i] = 0;
	ReadSeq(pfile_for, nseq, len, seq_real, olen, octa_av, octa1, octa, octa_size, octa_pro1, len_peak_max, outlog);
	for (i = 0; i < octa_size; i++)octa_rat[i] = log10(octa_av[i]);
	for (i = 0; i < octa_size; i++)octa_av[i] = 0;
	for (i = 0; i < octa_size; i++)octa1[i] = 0;
	for (i = 0; i < 16; i++)pexp[i] = 0;
	ReadSeqBack(pfile_back, nseqb, lenb, olen, reg_max, dav, dcv, octa_av, pexp, octa1, octa, octa_size, len_max, len_peak_max, outlog);
	for (i = 0; i < octa_size; i++)octa_rat[i] -= log10(octa_av[i]);
	/*fprintf(outlog,"Octa_rat");
	for (i = 0; i < octa_size; i++)
	{
		if (i % 16 == 0)fprintf(outlog,"\n%d", i);
		fprintf(outlog,"\t%g", octa_rat[i]);
	}
	fprintf(outlog,"\n");*/
	for (i = 0; i < nseq; i++)
	{
		int leni = len[i] - octa + 1;
		for (n = 0; n < leni; n++)
		{
			int i_sost = (int)octa_pro1[i][n];
			if (i_sost != -1)octa_pro1[i][n] = octa_rat[i_sost];
			else octa_pro1[i][n] = 0;
		}
	}
	//exit(1);
	for (j = 0; j < reg_max; j++)
	{
		fprintf(outlog, "AvR%d\t", j + 1);
		for (i = 0; i < 16; i++)fprintf(outlog, "%.4f\t", dav[j][i]);
		fprintf(outlog, "\n");
	}
	for (j = 0; j < reg_max; j++)
	{
		fprintf(outlog, "CovR%d\t", j + 1);
		for (i = 0; i < 16; i++)fprintf(outlog, "%.4f\t", dcv[j][i]);
		fprintf(outlog, "\n");
	}
	{
		char word[] = "acgt";
		GetWords(2, 0, 16, word);
	}
	strcpy(filef, file_for);
	strcat(filef, "_");
	strcpy(file_out, file_for);
	strcat(file_out, "_andy");
	for (i = 0; i < MEGE; i++)
	{
		pop[i].mem_in(nseq);
	}
	for (i = 0; i < 2; i++)
	{
		det2[i].mem_in(nseq);
	}
	det1.mem_in(nseq);
	qp = new double[nseq];
	if (qp == NULL) { puts("Out of memory..."); exit(1); }
	for (i = 0; i < nseq; i++)qp[i] = 0;
	thr_octa = new double[nseq];
	if (thr_octa == NULL) { puts("Out of memory..."); exit(1); }
	len_octa = new int[nseq];
	if (len_octa == NULL) { puts("Out of memory..."); exit(1); }
	octa_prowb = new int* [nseq];
	if (octa_prowb == NULL) return -1;
	for (i = 0; i < nseq; i++)
	{
		octa_prowb[i] = new int[len[i]];
		if (octa_prowb[i] == NULL) return -1;
	}
	octa_prows = new int* [nseq];
	if (octa_prows == NULL) return -1;
	for (i = 0; i < nseq; i++)
	{
		octa_prows[i] = new int[len[i]];
		if (octa_prows[i] == NULL) return -1;
	}
	for (i = 0; i < nseq; i++)for (n = 0; n < len[i]; n++)octa_prowb[i][n] = -1;
	int len_tot = 0, len_wei = 0;
	for (i = 0; i < nseq; i++)
	{
		int leni = len[i] - olen + 1;
		int odif = olen - octa + 1;
		for (n = 0; n < leni; n++)
		{
			octa_prow[i][n] = 0;
			for (k = 0; k < odif; k++)octa_prow[i][n] += octa_pro1[i][n + k];
			octa_prow[i][n] /= odif;
		}
		len_tot += leni;
		int half = leni / 3 - 1;
		{
			for (n = 0; n < leni; n++)octa_pro1p[n] = octa_prow[i][n];
			qsort(octa_pro1p, leni, sizeof(double), compare_qq2);
			thr_octa[i] = octa_pro1p[half];
		}
		double maxw = 0;
		k = 0;
		for (n = 0; n < leni; n++)
		{
			double dw = octa_prow[i][n] - thr_octa[i];
			if (dw >= 0)
			{
				octa_prowb[i][k] = n;
				k++;
				if (dw > maxw)dw = maxw;
			}
		}
		len_octa[i] = k;
		len_wei += len_octa[i];
		double koef;
		if (maxw > 0) { koef = Max(1, 5 / maxw); }
		else koef = 1;
		octa_prows[i][0] = 1 + (int)(koef * (octa_prow[i][octa_prowb[i][0]] - thr_octa[i]));
		for (n = 1; n < len_octa[i]; n++)
		{
			octa_prows[i][n] = octa_prows[i][n - 1] + 1 + (int)(koef * (octa_prow[i][octa_prowb[i][n]] - thr_octa[i]));
		}
	}
	int big_exit1 = 1;// local exit (separ +-) global exit (separation do not exceeded the previous run)
	double fit_prev, fit_after_mut;
	int cnt_count = 0;
	//Test(peak_real[0],len,0,2);		
	fprintf(outlog, "\n%s\tTrain\tFractHoxa %f\t", file_for, (double)len_wei / len_tot);
	fprintf(outlog, "Ndi %d\tDeg %d\tEli %d BE1 %d\n", size, MEGE, MEGE, big_exit1);
	//initiation		
	int gen = 0;
	int check_lpd = 1;
	int pair_all;
	int** pair_d; // pair 1 & 2 population ranks of individuals in each cell				
	int* pair_take;
	int pair_all_max = 0;
	double jwei0[3];
	{
		//= { 1.02, 1.05, 1.1 };
		int rwei[MEGE];
		double jwei1 = pow(double(20), double(2) / (MEGE + 2));//= jwei0[2];
		jwei0[2] = jwei1;
		jwei0[1] = 1 + (jwei1 - 1) / 2;
		jwei0[0] = 1 + (jwei1 - 1) / 5;
		int jmax = MEGE / 2 - 1, jmax1 = jmax - 1;
		rwei[jmax] = 2;
		for (j = jmax1; j >= 0; j--)
		{
			rwei[j] = (int)(rwei[jmax] * jwei1);
			jwei1 *= jwei0[2];
		}
		printf("RecWei\t");
		for (j = jmax; j >= 0; j--)printf("%d %d\t", j, rwei[j]);
		printf("\n");
		for (j = jmax; j >= 0; j--)
		{
			for (k = j; k < MEGE; k++)
			{
				for (m = 0; m < rwei[j]; m++)
				{
					pair_all_max += 2;
				}
			}
		}
	}
	fprintf(outlog, "Total pairs reserved %d\n", pair_all_max);
	pair_take = new int[pair_all_max];
	if (pair_take == NULL) return -1;
	pair_d = new int* [pair_all_max];
	if (pair_d == NULL) return -1;
	for (m = 0; m < pair_all_max; m++)
	{
		pair_d[m] = new int[2];
		if (pair_d[m] == NULL) return -1;
	}
	int success_r, success_m;
	int success_r1[NREC];
	//PARAMETERS SETTING																
	int stop_oi[MEGE], stop_li[MEGE], stop_pi[MEGE];
	for (i = 0; i < MEGE; i++)stop_oi[i] = stop_li[i] = stop_pi[i] = 0;
	int mege_h;
	//int restart_full = 0, restart_half = 0;
	do
	{
		int success_o, success_l, success_p;
		if (big_exit1 == 1)
		{
			for (i = 0; i < MEGE; i++)
			{
				int err = 1;
				//			fprintf(outlog,"Poo Ini1 %d\n", i);							
				do
				{
					err = 1;
					int gom = 0;
					det1.init_rand_hoxa(nseq, olen, size, reg_max, len_octa, len, octa_prowb, octa_prows);
					if (det1.check(0, reg_max, outlog) == -1)
					{
						det1.check(0, reg_max, outlog);
						fprintf(outlog, "Population error!\n");
						det1.print_all(reg_max, nseq, outlog);
						exit(1);
					}
					else
					{
						for (m = 0; m < i; m++)
						{
							gom = GomTown(det1, pop[m], nseq, 1);
							if (gom == -1)
							{
								break;
							}
						}
					}
					if (gom != -1)
					{
						//							fprintf(outlog,"Poo Ini2 %d\n", i);
						EvalMahFIT(&det1, nseq, octa, infc, seq_real, olen, dav, dcv, octa_prow, pexp, len, outlog);//hoxa_wei,qps, frp, octa_rat,
						det1.print_sta(reg_max, outlog);
						det1.get_copy(&pop[i], nseq, reg_max);
						//						fprintf(outlog,"Poo Ini3 %d\n", i);
						err = 0;
					}
				} while (err == 1);
			}
			big_exit1 = 0;
			qsort((void*)(&pop[0]), MEGE, sizeof(pop[0]), compare_pop);
			//pop[iter][0].print_all(reg_max,nseq);
		}
		/*for (i = 0; i<MEGE; i++)
		{
		//fprintf(outlog,"After %d\n",i+1);
		//	pop[i].print_all(reg_max,nseq);
		if (pop[i].check(0, reg_max) == -1)
		{
		pop[i].print_all(reg_max, nseq);
		fprintf(outlog,"Population error!\n");
		exit(1);
		}
		}*/
		fit_prev = pop[0].fit;
		//pop[0].print_all(reg_max,nseq);													
		success_o = success_l = success_p = success_m = 0;
		double ratio2_gen0 = 0.01, ratio_rec_cycle, ratio_mut_cycle = 0.0001;
		int step, step_max, step_max_tot = 0;
		int elit_rec;
		int knseq = Max(100, nseq);
		int sr_step = 400 * knseq;//4000000 if nseq = 1000		;
		int n_rec_cycle_max = 4000 * knseq;//4000000 if nseq = 1000		
		double jwei;
		if (gen == 0)
		{
			jwei = jwei0[0];
			step = 2000;
			step_max = 20 * knseq; //20000 if nseq = 1000			
			elit_rec = MEGE / 5;			
			mege_h = MEGE;
			ratio_rec_cycle = 0.0025;
			sr_step = 400 * knseq;
		}
		else
		{
			if (gen == 1)
			{
				jwei = jwei0[1];
				step = 4000;
				step_max = 200 * knseq;//500000 if nseq = 500				
				elit_rec = MEGE / 5;
				mege_h = MEGE;
				ratio_rec_cycle = 0.001;
				sr_step = 1000 * knseq;
			}
			else
			{
				jwei = jwei0[2];
				step = 5000;
				step_max = 200 * knseq;//100000 if nseq = 500				
				elit_rec = MEGE/5;				
				mege_h = MEGE;
				ratio_rec_cycle = 0.001;
				sr_step = 1000 * knseq;
			}
		}
		if (step_max < step)step_max = step;
		double ratio_thr = 1 / (double)step;
		double ratio_thr_r[2];
		for (i = 0; i < 2; i++)ratio_thr_r[i] = ratio_thr;
		int step2 = 2 * step;
		int mege_h1 = mege_h - 1;
		//mutations						
		int n_mut_tot = 0, success_m_tot = 0, mdo = 1;
		int asuccess[NMUT], atry[NMUT];
		for (k = 0; k < NMUT; k++)asuccess[k] = atry[k] = 0;
		int success_mi[MEGE], try_mi[MEGE];
		int success_ri[MEGE];
		if (gen == 0)for (i = 0; i < mege_h; i++)success_ri[i] = 1;
		double exp_rec_rate[MEGE];
		for (i = 0; i < mege_h; i++)exp_rec_rate[i] = 0;
		for (i = 0; i < mege_h; i++)success_mi[i] = try_mi[i] = 0;
		double sm0_rate = 0.001;
		int mut_jump = 0;
		do
		{
			step_max_tot += step_max;
			for (i = 0; i < mege_h; i++)pop[i].sm = pop[i].tm = 0;
			if (gen == 0)fprintf(outlog, "Mut cycle %d\n", step_max_tot / step_max);
			for (i = 0; i < mege_h; i++)
			{
				if (stop_pi[i] == 0 || (stop_li[i] == 0 || stop_oi[i] == 0))
				{													
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
					double fit_mut_prev0 = pop[i].fit;
					int n_mut_min_cyc = step;
					int m_iter = 0;
					while (ratio[2] > ratio_thr || (ratio[0] > ratio_thr || ratio[1] > ratio_thr))
					{
						pop[i].get_copy(&det1, nseq, reg_max);
						int sm;
						/*{
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
						}*/
						if (stop_oi[i] == 1)
						{
							if (stop_li[i] == 1)sm = 2;
							else
							{
								sm = rand() % 2;// 1 or 2
								sm++;
							}
						}
						else
						{
							if (stop_li[i] == 1)
							{
								sm = rand() % 2;
								if (sm == 1)sm++;// 0 or 2
							}
							else sm = rand() % 3;
						}
						int muto, muto1;
						int npeak = -1, nori = -1, npos = -1;
						if (sm == 0)muto = MutOlig0(&det1);
						else
						{
							if (sm == 1)muto = MutCry0(&det1, olen, reg_max);
							else
							{
								muto = MutRegShiftHoxaW(&det1, nseq, npeak, nori, npos, len_octa, octa_prowb, octa_prows);
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
									gom = GomTown(det1, pop[m], nseq, check_lpd);
									if (gom == -1) { break; }
								}
							}
						}
						//if (det1.check(0, reg_max) == -1){ fprintf(outlog,"Population error!\n"); exit(1); }
						//det.print_all();										
						double ret = 0, dd = 0;
						if (gom != -1)
						{
							n_mut_here++;
							ret = EvalMahFIT(&det1, nseq, octa, infc, seq_real, olen, dav, dcv, octa_prow, pexp, len, outlog);
							if (ret > 0)
							{
								step_try[sm]++;
								try_mi[i]++;
								dd = det1.fit / pop[i].fit;
								pop[i].tm++;
								if (dd > 1)
								{
									pop[i].sm++;
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
							}
						}
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
							int step_success_tot = step_success[0] + step_success[1] + step_success[2];
							if (step_success_tot < n_mut_min_cyc)n_mut_min_cyc = step_success_tot;
							double ratio_per_cycle = pop[i].fit / fit_mut_prev0 - 1;
							if (ratio[0] <= ratio_thr)stop_oi[i] = 1;
							if (ratio[1] <= ratio_thr)stop_li[i] = 1;
							if (ratio[2] <= ratio_thr)stop_pi[i] = 1;
							for (k = 0; k < 3; k++)atry[k] += step_try[k];
							for (k = 0; k < 3; k++)asuccess[k] += step_success[k];
							if (try_mi[i] != 0)exp_rec_rate[i] = (double)success_mi[i] / try_mi[i];
							else exp_rec_rate[i] = 0;
							//	if (gen >= 1)fprintf(outlog, "Step%d M%d %d,%d,%d Try %d,%d,%d Min %d M %f H %g Fit %f Ratio %f RatioSco %f\n", m_iter, i + 1, step_success[0], step_success[1], step_success[2], step_try[0], step_try[1], step_try[2], n_mut_min_cyc, pop[i].mah, pop[i].fpr, pop[i].fit, ratio[2], ratio_per_cycle);
							success_mi[i] = try_mi[i] = 0;
							for (k = 0; k < 3; k++)step_try[k] = step_success[k] = 0;
							if (n_mut_here >= step_max)
							{
								if (gen == 0)break;
								else if (stop_li[i] == 1)break;
							}
							if (ratio_per_cycle <= ratio_mut_cycle)
							{
								//	fprintf(outlog,"Too small score growth %f\n", ratio_per_cycle);
								stop_pi[i] = 1;
								break;
							}
							fit_mut_prev0 = pop[i].fit;
						}
						//fprintf(outlog,"\n");
					}					
					success_o += success_o_local;
					success_l += success_l_local;
					success_p += success_p_local;
					success_m_tot += success_m_local;
				//	if ((gen == 0 && i == 0) || gen > 0)
					{
						fprintf(outlog, "M %d %d,%d,%d = %d Try %d M %f H %g F %f", i + 1, success_o_local, success_l_local, success_p_local, success_m_local, n_mut_here, pop[i].mah, pop[i].fpr, pop[i].fit);
						fprintf(outlog, "\n");
					}
					/*if (gen == 0 && i == 0)
					{
						fprintf(outlog, "\tOLP %d%d%d RatioP %f", stop_oi[i], stop_li[i], stop_pi[i], ratio[2]);
					}*/
					n_mut_tot += n_mut_here;
					if (success_m_local != 0)success_m++;
				}
			}
			mut_jump = 0;
			for (i = 1; i < mege_h; i++)
			{
				if (pop[i].fit >= pop[0].fit)mut_jump++;
			}
			qsort((void*)(&pop[0]), mege_h, sizeof(pop[0]), compare_pop);			
			{
				int trym = 0, sucm = 0;
				for (i = 0; i < mege_h; i++)
				{
					sucm += pop[i].sm;
					trym += pop[i].tm;
				}
				sm0_rate = (double)sucm / trym;
			}
			if (gen > 0)
			{
				mdo = 0;
				break;
			}
			else
			{
				mdo = 0;
				int mege_half = mege_h / 2;
				if (mut_jump >= mege_half)
				{
					mdo = 1;
					fprintf(outlog, "Mut jump %d\tSM rate %g\n", mut_jump, sm0_rate);
				}
				else
				{
					int m_nostoplo = 0, m_stopp = 0;
					for (i = 0; i < mege_h; i++)
					{
						if ((stop_li[i] == 0 || stop_oi[i] == 0) && stop_pi[i] == 0)
						{
							m_nostoplo++;
						}
						if (stop_pi[i] == 1)
						{
							m_stopp++;
						}
					}
					double rat0 = 0;
					if (atry[2] > 0)
					{
						rat0 = (double)asuccess[2] / atry[2];
					}
					else mdo = 0;
					if (m_nostoplo > mege_half || rat0 > ratio2_gen0)mdo = 1;
					if (m_stopp >= mege_half)mdo = 0;
					fprintf(outlog, "Common RatioP %f\tMut jump %d\tyLOnP %d yP %d\tSM rate %g\t", rat0, mut_jump, m_nostoplo, m_stopp, sm0_rate);
					tnow = time(NULL);
					fprintf(outlog, "%s", ctime(&tnow));
				}
			}
		} while (mdo == 1);
		qsort((void*)(&pop[0]), mege_h, sizeof(pop[0]), compare_pop);
		if (gen == 0)success_m /= (step_max_tot / step_max);
		if (gen <= 1)ratio_thr_r[0] = ratio_thr;		//OL			
		else ratio_thr_r[0] = 1;		
		for (i = 0; i < mege_h; i++)pop[i].sr = pop[i].tr = 0;
		{
			int recount = 0;
			double resum = 0;
			for (i = 0; i < mege_h - 1; i++)
			{
				double e1 = Max(ratio_thr, exp_rec_rate[i]);
				double e2 = Max(ratio_thr, exp_rec_rate[i + 1]);
				resum += success_ri[i] * (e1 + e2) / 2;
				recount += success_ri[i];
			}
			if (recount != 0)ratio_thr_r[1] = resum / recount;
			else ratio_thr_r[1] = ratio_thr;
		}
		//recombinations											
		fprintf(outlog, "Total mut %d\tSM rate %g\t Mut jump %d\n", n_mut_tot, sm0_rate, mut_jump);
		int loc_rec, loc_rec_tot = 0;
		fit_after_mut = pop[0].fit;
		for (m = 0; m < NREC; m++)success_r1[m] = 0;
		int rwei[MEGE];
		{
			double jwei1 = jwei;
			int jmax = mege_h / 2 - 1, jmax1 = jmax - 1;
			rwei[jmax] = 2;
			for (j = jmax1; j >= 0; j--)
			{
				rwei[j] = (int)(rwei[jmax] * jwei1);
				jwei1 *= jwei;
			}
			for (j = jmax; j >= 0; j--)printf("%d %d\t", j, rwei[j]);
			printf("\n");
			pair_all = 0;
			{
				for (j = jmax; j >= 0; j--)
				{
					for (k = j + 1; k < mege_h; k++)
					{
						for (m = 0; m < rwei[j]; m++)
						{
							pair_d[pair_all][0] = j;
							pair_d[pair_all][1] = k;
							pair_all++;
							pair_d[pair_all][0] = k;
							pair_d[pair_all][1] = j;
							pair_all++;
						}
					}
				}
			}
			fprintf(outlog, "\nTotal pairs %d Gen %d\n", pair_all, gen + 1);
		}
		for (k = 0; k < pair_all; k++)pair_take[k] = k;
		int n_rec_cycle, sr_step_cycle = sr_step;
		if (gen <= 1) { n_rec_cycle = Min(n_mut_tot / 4, n_rec_cycle_max); }
		else n_rec_cycle = n_rec_cycle_max;
		n_rec_cycle /= pair_all;
		if (n_rec_cycle < 1)n_rec_cycle = 1;
		sr_step_cycle /= pair_all;
		int rec_jump[MEGE];
		for (k = 0; k < mege_h; k++)rec_jump[k] = 0;
		double fit_pre_rec[MEGE];
		for (k = 0; k < mege_h; k++)fit_pre_rec[k] = pop[k].fit;
		int sr;
		double ratio_r[2] = { 1,1 };
		int step_rtry[NREC], step_rsuccess[NREC];
		for (k = 0; k < NREC; k++)step_rtry[k] = step_rsuccess[k] = 0;
		int success_r_cycle = 0;
		fprintf(outlog, "Rec %d cycles of %d tries\tTotal %d\tRatioThrOL %.5f RatioThrP %.5f StepR %d\t", n_rec_cycle, pair_all, n_rec_cycle * pair_all, ratio_thr_r[0], ratio_thr_r[1], sr_step_cycle);
		fprintf(outlog, "RE %d", elit_rec);
		fprintf(outlog, "\n");
		double fit_rec_prev0 = pop[0].fit;
		for (i = 0; i < mege_h; i++)success_ri[i] = 0;
		loc_rec = 0;
		for (i = 0; i < mege_h; i++)pop[i].sr = pop[i].tr = 0;
		double fit_rec_prev = pop[elit_rec].fit;
		double mah_rec_prev = pop[0].mah;
		double fpr_rec_prev = pop[0].fpr;
		sr = 0;
		//for (sr = 1; sr <= n_rec_cycle; sr++)
		do
		{
			sr++;
			BigMixI(pair_take, pair_all);
			for (k = 0; k < pair_all; k++)
			{
				int kk[2];
				for (m = 0; m < 2; m++)kk[m] = pair_d[pair_take[k]][m];
				for (m = 0; m < 2; m++)
				{
					pop[kk[m]].get_copy(&det2[m], nseq, reg_max);
				}
				double fit_parent_max;
				fit_parent_max = Max(det2[0].fit, det2[1].fit);
				int nsn[2], nsi[2], num[2];
				int sumw = 0, r_cy;
				int kk01 = pair_take[k];
				r_cy = rand() % 6;
				switch (r_cy)
				{
				case(0):
				{
					if (Reco2_Original(&det2[0], &det2[1], nsi, num) == -1)continue;
					else break;
				}
				case(1):
				{
					if (Reco2_Economic(&det2[0], &det2[1], nsi, num) == -1)continue;
					else break;
				}
				case(2):
				{
					if (Reco2_One_dinucleotide_local(&det2[0], &det2[1], olen, reg_max) == -1)continue;
					else break;
				}
				case(3):
				{
					if (Reco2_One_dinucleotide_full(&det2[0], &det2[1]) == -1)continue;
					else break;
				}
				case(4):
				{
					if (Reco2_Two_dinucleotides_full(&det2[0], &det2[1], nseq, reg_max) == -1)continue;
					else break;
				}
				case(5):
				{
					if (Reco2Peak(&det2[0], &det2[1], nseq) == -1)continue;
					else break;
				}
				}
				//fprintf(outlog,"out Rec %d\n",r_cy);
				if (r_cy <= 1)
				{
					for (m = 0; m < 2; m++)
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
				fprintf(outlog,"Population error!\n");
				exit(1);
				}
				}*/
				double dd_r[2] = { 0, 0 };
				int gom[2] = { 0, 0 };
				{
					for (m = 0; m < 2; m++)
					{
						for (int t = 0; t < mege_h; t++)
						{
							if (t == pair_d[k][m])continue;
							gom[m] = GomTown(det2[m], pop[t], nseq, check_lpd);
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
						det2[m].fit = EvalMahFIT(&det2[m], nseq, octa, infc, seq_real, olen, dav, dcv, octa_prow, pexp, len, outlog);
					}
					step_rtry[r_cy]++;
					for (m = 0; m < 2; m++)pop[kk[m]].tr++;
					{
						for (m = 0; m < 2; m++)
						{
							if (det2[m].fit > pop[kk[m]].fit)
							{
								pop[kk[m]].sr++;
								success_r1[r_cy]++;
								step_rsuccess[r_cy]++;
								loc_rec++;
								rec_jump[kk[m]]++;
								success_ri[kk[m]]++;
								success_r_cycle++;
								det2[m].get_copy(&pop[kk[m]], nseq, reg_max);
							}
						}
					}
				}
			}
			success_r = 0;
			for (i = 0; i < mege_h; i++)if (success_ri[i] > 0)success_r++;
			if (sr % sr_step_cycle == 0)
			{
				int stepr;
				stepr = step_rtry[0] + step_rtry[1] + step_rtry[2] + step_rtry[3] + step_rtry[4];
				if (stepr > 0)ratio_r[0] = (double)(step_rsuccess[0] + step_rsuccess[1] + step_rsuccess[2] + step_rsuccess[3] + step_rsuccess[4]) / stepr;
				else ratio_r[0] = 0;
				stepr = step_rtry[5];
				if (stepr > 0)ratio_r[1] = (double)(step_rsuccess[5]) / stepr;
				else ratio_r[1] = 0;
				qsort((void*)(&pop[0]), mege_h, sizeof(pop[0]), compare_pop);
				int tryr = 0, sucr = 0, rec_jump1 = 0;
				for (i = 0; i < mege_h; i++)
				{
					if (rec_jump[i] > 0)
					{
						rec_jump1++;
						sucr += pop[i].sr;
						tryr += pop[i].tr;
					}
					else break;
				}
				double sr0_rate;
				if (tryr > 0) sr0_rate = (double)sucr / tryr;
				else sr0_rate = 0;
				fit_rec_prev = pop[elit_rec].fit;
				loc_rec_tot += loc_rec;
				fprintf(outlog, "Rec %d: %d %d,%d,%d,%d,%d,%d %d OL %f P %f M %f FP %g F %f ", sr * pair_all, success_r, success_r1[0], success_r1[1], success_r1[2], success_r1[3], success_r1[4], success_r1[5], success_r_cycle, ratio_r[0], ratio_r[1], pop[0].mah, pop[0].fpr, pop[0].fit);
				loc_rec_tot += loc_rec;
				double ratio_per_cycle = pop[0].fit / fit_rec_prev0 - 1;
				fprintf(outlog, "L%d RE %f SR rate %f Rec Jump %d RatioSco %f\n", loc_rec_tot, fit_rec_prev, sr0_rate, rec_jump1, ratio_per_cycle);
				for (k = 0; k <= elit_rec; k++)fprintf(outlog, "%d M %f H %f F %f\n", k + 1, pop[k].mah, pop[k].fpr, pop[k].fit);
				for (k = 0; k < mege_h; k++)fprintf(outlog, "%d ", rec_jump[k]);
				fprintf(outlog, "\n");
				tnow = time(NULL);
				fprintf(outlog, "%s", ctime(&tnow));
				for (k = 0; k < mege_h; k++)rec_jump[k] = 0;
				if (loc_rec == 0)break;
				if (ratio_r[0] < ratio_thr_r[0] && ratio_r[1] < ratio_thr_r[1])break;
				if (sr0_rate <= sm0_rate)break;
				if (ratio_per_cycle <= ratio_rec_cycle)break;
				if (sr > n_rec_cycle && rec_jump1 <= elit_rec + 1)break;
				for (m = 0; m < NREC; m++)step_rtry[m] = step_rsuccess[m] = 0;
				success_r_cycle = 0;
				loc_rec = 0;
				fit_rec_prev0 = pop[0].fit;
			}
		} while (loc_rec >= 0);
		qsort((void*)(&pop[0]), mege_h, sizeof(pop[0]), compare_pop);
		for (k = 1; k < mege_h; k++)
		{
			if (stop_pi[k] == 1 && fit_pre_rec[k] > pop[k].fit)stop_pi[k] = 0;
		}
		double change_level = pop[0].fit / fit_prev - 1;
		double change_level_rec = pop[0].fit / fit_after_mut - 1;
		double change_level_mut = fit_after_mut / fit_prev - 1;
		double exit_1st = 0.05, exit_2nd = 0.01;
		if (change_level < exit_2nd && gen > 1)// 3rd iteration at least
		{
			printf("Go out!\n");
			big_exit1 = 1;
		}
		else
		{
			fprintf(outlog, "Continue\n");
			if (gen > 1 && change_level < exit_2nd)
			{
				/*int best_alive = (int)(mege_h / (3 + 2 * restart_half));//for mege = 30 -> 10,6,4,3,2,2,2
				if (best_alive < 2)best_alive = 2;
				int n_peaks = 20;
				restart_full++;
				fprintf(outlog, "Restart Full %d\n", restart_full);
				for (i = 0; i < mege_h; i++)stop_oi[i] = stop_li[i] = 0;
				for (i = best_alive; i < mege_h; i++)
				{
					int best_num = i % best_alive;
					//	if (GomTown2(pop[iter][j][best_num], pop[iter][j][i]) == -1)continue;
					pop[best_num].get_copy(&pop[i], nseq, reg_max);
					pop[i].init_rand_part_hoxa(nseq, n_peaks, olen, len_octa, len, octa_prowb, octa_prows);
					EvalMahFIT(&pop[i], nseq, octa, seq_real, olen, dav, dcv, octa_prow, len, outlog);//qps, frp, octa_rat,
				}*/
				qsort((void*)(&pop[0]), mege_h, sizeof(pop[0]), compare_pop);
				for (j = 0; j < mege_h; j++)stop_pi[i] = stop_oi[i] = stop_li[i] = 0;
			}
			/*if (restart_full == 0)
			{
				int best_alive = (int)(mege_h / (3 + 2 * restart_half));//for mege = 30 -> 10,6,4,3,2,2,2
				if (best_alive < 2)best_alive = 2;
				restart_half++;
				fprintf(outlog, "Restart Half %d\n", restart_half);
				int n_peaks = 20;
				for (i = best_alive; i < mege_h; i++)
				{
					int best_num = i % best_alive;
					pop[best_num].get_copy(&pop[i], nseq, reg_max);
					pop[i].init_rand_part_hoxa(nseq, n_peaks, olen, len_octa, len, octa_prowb, octa_prows);
					EvalMahFIT(&pop[i], nseq, octa, seq_real, olen, dav, dcv, octa_prow, len, outlog);
				}
				qsort((void*)(&pop[0]), mege_h, sizeof(pop[0]), compare_pop);
				for (j = 0; j < mege_h; j++)stop_pi[i] = 0;
				if (restart_half == 1)
				{
					for (j = 0; j < mege_h; j++)stop_oi[i] = stop_li[i] = 0;
				}
				else
				{
					for (j = best_alive; j < mege_h; j++)stop_oi[i] = stop_li[i] = 0;
				}
			}*/
		}
		gen++;
		fit_prev = pop[0].fit;
		//if (change_level<GA_EXIT)gen1++;
		fprintf(outlog, "Gen %d Fit %.5f Rat %.5f RatM %.5f RatR %.5f ", gen, pop[0].fit, change_level, change_level_mut, change_level_rec);
		//	if(restart==1)fprintf(outlog,"W %d WP %d AvS %.1f AvP %.1f", wei_max, wei_max_pos, wei_peak_av, wei_pos_av);
		fprintf(outlog, "M %d %d,%d,%d R %d ", success_m, success_o, success_l, success_p, success_r1[4]);
		//	if (restart == 0)fprintf(outlog,"LM %d LR %d ", loc_mut_tot, loc_rec_tot);
		{
			//	fprintf(outlog,"R %d %d,%d,%d,%d,%d ", success_r, success_r1[0], success_r1[1], success_r1[2], success_r1[3], success_r1[4]);
			//	fprintf(outlog,"%d,%d,%d,%d,%d ", recw[0][0], recw[1][0], recw[2][0], recw[3][0], recw[4][0]);
			//	fprintf(outlog,"%d,%d,%d,%d,%d ", recw[0][pair_all1], recw[1][pair_all1], recw[2][pair_all1], recw[3][pair_all1], recw[4][pair_all1]);				
		}
		{
			int sumr = 0;
			for (i = 0; i < mege_h; i++)sumr += success_ri[i];
			if (sumr > 0)for (i = 0; i < mege_h; i++)fprintf(outlog, "%.2f ", 100 * (double)success_ri[i] / sumr);
		}
		tnow = time(NULL);
		fprintf(outlog, "%s", ctime(&tnow));
		fit_prev = pop[0].fit;
		//	if (restart == 0)for (i = 0; i < 2; i++)pop[i].print_all(reg_max, nseq);
		//if (gen == 1)for (i = 0; i < mege_h; i++)stop_li[i] = stop_oi[i] = 1;
		for (i = 0; i < mege_h; i++)stop_pi[i] = 0;
	} while (big_exit1 == 0);
	char name[500];
	for (i = 0;; i++)
	{
		if (file_for[i] == '.') { name[i] = '\0'; break; }
		if (file_for[i] == '\n') { name[i] = '\0'; break; }
		if (file_for[i] == '\0') { name[i] = '\0'; break; }
		name[i] = file_for[i];
	}
	char extmat0[] = "_mat";
	char extmap0[] = "_loc";
	char extmat[6], extmap[6];
	char file_base[500];
	strcpy(file_base, path_out);
	strcat(file_base, name);
	for (i = 0; i < MEGE; i++)
	{
		char bufi[5];
		memset(bufi, '\0', sizeof(bufi));
		sprintf(bufi, "%d", i + 1);
		memset(extmat, '\0', sizeof(extmat));
		memset(extmap, '\0', sizeof(extmap));
		strcpy(extmat, extmat0);
		strcpy(extmap, extmap0);
		strcat(extmat, bufi);
		strcat(extmap, bufi);
		EvalMahFITTrain(&pop[i], nseq, seq_real, file_for, olen, reg_max, len, dav, dcv, outlog, file_base, extmap, extmat, name);
	}
	fprintf(outlog, "Go out big cycle ");
	big_exit1 = 1;
	fclose(outlog);
	for (i = 0; i < MEGE; i++)pop[i].mem_out();
	for (i = 0; i < 2; i++)det2[i].mem_out();
	det1.mem_out();
	for (k = 0; k < 2; k++)
	{
		for (i = 0; i < nseq; i++)
		{
			delete[] seq_real[k][i];
		}
		delete[] seq_real[k];
	}
	delete[] seq_real;
	for (i = 0; i < pair_all_max; i++)
	{
		delete[] pair_d[i];
	}
	delete[] pair_d;	
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
	delete[] qp;
	delete[] octa_rat;
	delete[] octa_av;
	delete[] octa1;
	for (i = 0; i < nseq; i++)
	{
		delete[] octa_pro1[i];
	}
	delete[] octa_pro1;

	for (i = 0; i < nseq; i++)
	{
		delete[] octa_prow[i];
	}
	delete[] octa_prow;
	for (i = 0; i < nseq; i++)
	{
		delete[] octa_prowb[i];
	}
	delete[] octa_prowb;
	delete[] octa_pro1p;
	return 0;
}

