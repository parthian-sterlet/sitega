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
#define SEQLEN 5000
#define LPDLEN 6 //max LPD length
#define MOTLEN 20 //max motif len
#define CELL 4//no. of cell populations
#define MEGE 20//population size
#define NMUT 3
#define NREC 6
#define POPSIZE 80
//#define CENT 100
//#define GA_EXIT 0.01

double uw[POPSIZE][POPSIZE];

struct ss {
	int num;
	char oli[3];
}s[16];

int compare_qq(const void* X1, const void* X2)//increase
{
	double X = (*(double*)X1 - *(double*)X2);
	if (X > 0)return 1;
	if (X < 0)return -1;
	return 0;
}
int compare_qq2(const void* X1, const void* X2)//decrease
{
	double X = (*(double*)X1 - *(double*)X2);
	if (X > 0)return -1;
	if (X < 0)return 1;
	return 0;
}
struct uno {
	int sta;
	int end;
	int num;
	void get_copy(uno* a);
	void print_all(FILE* out);
	void print_all2();
};
void uno::get_copy(uno* a)
{
	a->num = num;
	a->sta = sta;
	a->end = end;
};
void uno::print_all(FILE* out)
{
	fprintf(out, "[%d;%d]%s\t", sta, end, s[num].oli);
}
void uno::print_all2()
{
	printf("[%d;%d]%s\t", sta, end, s[num].oli);
}
struct town {
	uno tot[POPSIZE];
	int deg[16];
	int size;
	int* pos;// pozicii na4al okon
	//int *posc;// pozicii na4al okon
	int* ori;// DNA strand 0,1
	double fit;
	double mah;
	double inf;
	double fpr;
	double aucroc;
	double auprc;
	int sm;// success mut
	int sr;// success rec
	int tm;// try mut
	int tr;// try rec
	int odg[LPDLEN + 1];
	void get_copy(town* a, int nseq, int reg_max);
	//void init_rand(int nseq, int *len, int oln, int rsize, int reg_max);
	void init_rand_hoxa(int nseq, int oln, int rsize, int reg_max, int* len_octa, int* len, int** octa_prowb, int** octa_prows);
	//void init_rand_part(int nseq, int *len, int oln, int *xporti, int nind);
	void init_rand_part_hoxa(int nseq, int* xporti, int nind, int olen, int* len_octa, int* len, int** octa_prowb, int** octa_prows);
	int init_add(uno last);
	int init_add_split(void);
	void init_zero(int olen);
	void Mix(int* a, int* b);
	void BigMix(int* d, int ln);
	void swap(int n1, int n2);
	int order(int n);
	void print_all(int reg_max, int nseq, FILE* out);
	void print_all2(int reg_max);
	void print_sta(int reg_max, FILE* out);
	void sort_all(void);
	int sum(int j);
	void fprint_all(char* file, char* add);
	void fprint_allfi_mat(char *file, int len, double sga_min, double sga_raz, double* buf);
	void fprint_seq(char* file, char* add, int olen, int nseq, char*** seq, double* best_sco, int* len, int* xportj, int cnt_shift);
	int check(int min, int max, FILE* out);
	int mem_in(int nseq);
	void mem_out(void);
}***pop, det1, det2[2];
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
	fpr = 1;
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
void town::fprint_seq(char* file, char *add, int olen, int nseq, char*** seq, double* best_sco, int* len, int *xportj, int cnt_shift)
{
	int i, j, k, x1, x2, m;
	char d[POPSIZE], dir[] = "+-";
	char file_out[500];
	strcpy(file_out, file);
	strcat(file_out, add);

	FILE* out;
	if ((out = fopen(file_out, "wt")) == NULL)
	{
		printf("Ouput file can't be opened!\n");
		exit(1);
	}
	for (m = 0; m < nseq; m++)
	{
		i = xportj[m];
		x1 = pos[i], x2 = pos[i] + olen;
		int cep = ori[i];
		k = 0;
		for (j = x1; j < x2; j++)d[k++] = seq[cep][i][j];
		d[olen] = '\0';
		if (cep == 1)
		{
			int x10 = x1, x20 = x2;
			x1 = len[i] - x20 + 1;
			x2 = len[i] - x10;
		}
		else
		{
			x1++;
		}
		fprintf(out, "%d\t%d\t%d\t%c\t%f\t%s\n", i + 1, x1, x2, dir[ori[i]], best_sco[cnt_shift+m], d);
	}
	fclose(out);
}
void town::fprint_allfi_mat(char* file, int len, double sga_min, double sga_raz, double* buf)
{
	int i;
	FILE* out;
	if ((out = fopen(file, "wt")) == NULL)
	{
		printf("Ouput file can't be opened!\n");
		exit(1);
	}
	fprintf(out, "Bootatrap\n");
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
void town::print_all(int reg_max, int nseq, FILE* out)
{
	int i;
	char strand[] = "+-";
	fprintf(out, "Mah %f FP %f Inf %f Fit %f S %d\t", mah, fpr, inf, fit, size);
	for (i = 0; i < 16; i++)fprintf(out, "%d ", deg[i]); fprintf(out, "\tLEN ");
	for (i = 0; i < reg_max; i++)fprintf(out, "%d ", odg[i]); fprintf(out, "\t");
	int size1 = size - 1;
	fprintf(out, " ");
	for (i = 0; i < size; i++)
	{
		tot[i].print_all(out);
		if (i == size1)fprintf(out, "\n");
		// else if(tot[i + 1].num != tot[i].num)fprintf(out,"\n");
		//	 fprintf(out," ");
	}
	// fprintf(out,"\n");

	/*for (i = 0; i < nseq; i++)
	{
		fprintf(out, "%2d", pos[i]);
		fprintf(out, "%c", strand[ori[i]]);
		fprintf(out, " ");
	}
	fprintf(out, "\n");*/
}
void town::print_all2(int reg_max)
{
	int i;
	char strand[] = "+-";
	printf("M %f FP %f F %f S %d\t", mah, fpr, fit, size);
	for (i = 0; i < 16; i++)printf("%d ", deg[i]); printf("\tLEN ");
	for (i = 0; i < reg_max; i++)printf("%d ", odg[i]); printf("\t");
	int size1 = size - 1;
	printf(" ");
	for (i = 0; i < size; i++)
	{
		tot[i].print_all2();
		if (i == size1)printf("\n");
		// else if(tot[i + 1].num != tot[i].num)fprintf(out,"\n");
		//	 fprintf(out," ");
	}
	// fprintf(out,"\n");

	/*for (i = 0; i < nseq; i++)
	{
		fprintf(out, "%2d", pos[i]);
		fprintf(out, "%c", strand[ori[i]]);
		fprintf(out, " ");
	}
	fprintf(out, "\n");*/
}
void town::print_sta(int reg_max, FILE* out)
{
	int i;
	fprintf(out, "M %f H %f I %f F %f\t", mah, fpr, inf, fit);
	for (i = 0; i < 16; i++)fprintf(out, "%2d", deg[i]); fprintf(out, "\t");
	for (i = 0; i < reg_max; i++)fprintf(out, "%2d ", odg[i]); fprintf(out, "\n");
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
	/*fprintf(out,"\n");
	fprintf(out,"old %d\t",n);
	tot[n].print_all();fprintf(out,"\n");
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
	//fprintf(out,"DIRECT %d\t",di);
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
	/*fprintf(out,"\n");
	fprintf(out,"new %d\t",ret);
	tot[ret].print_all();fprintf(out,"\n");
	print_all();
	fprintf(out,"\n");*/
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
			//		fprintf(out,"DO\t");
			//	for(j=0;j<oln1;j++)fprintf(out,"%d\t",take_pos[j]);fprintf(out,"\n");
			BigMixI(take_pos, oln1);
			//	fprintf(out,"PO\t");
			//	for(j=0;j<oln1;j++)fprintf(out,"%d\t",take_pos[j]);fprintf(out,"\n");
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
void town::init_rand_part_hoxa(int nseq, int* xporti, int nind, int olen, int* len_octa, int* len, int** octa_prowb, int** octa_prows)
{
	int i, j, k;

	fit = 0;
	fpr = 1;
	//int oln1 = olen - 1;
	for (j = 0; j < nseq; j++)
	{
		i = xporti[j];
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
int town::mem_in(int nseq)
{
	pos = new int[nseq];
	if (pos == NULL) return -1;
	//	posc = new int[nseq];
		//if (posc == NULL) return -1;
	ori = new int[nseq];
	if (ori == NULL) return -1;
	return 1;
}
void town::mem_out(void)
{
	delete[] pos;
	//	delete[] posc;
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
	char file_out[40];
	strcpy(file_out, file);
	strcat(file_out, add);
	if ((out = fopen(file_out, "at")) == NULL)
	{
		fprintf(out, "Ouput file can't be opened!\n");
		exit(1);
	}
	fprintf(out, "\t//%s\n", file);
	fprintf(out, "%d\t%.12f\t", size, fit);
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
		//fprintf(out,"%d\t%s\n",i,s[i].oli);
	}
}
void DelHole(char* str)
{
	char* hole;
	hole = strstr(str, "\n");
	if (hole != NULL) *hole = '\0';
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
int CheckStr(char* file, char* d, int n, int print, FILE* out)
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
		if (print == 1)fprintf(out, "File %s; sequence %d position %d (%c) bad. Sequence deleted!\n", file, n, i + 1, d[i]);
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
int BackMat(int size, FILE* out)
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
				fprintf(out, "\n1back %g(%d,%d) bij%g(%d,%d)\t", uw[i][j], i, j, b[i][j], i, j);
				//fprintf(out,"%g(%d)\n",uw[j][j],j);
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
				fprintf(out, "\n2back %g(%d,%d) bij%g(%d,%d)\t", uw[i][j], i, j, b[i][j], i, j);
				//fprintf(out,"%g(%d)\n",uw[j][j],j);
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
			fprintf(out, "\nback3 %g %d\t", uw[i][i], i);
			return -1;
		}
		for (j = 0; j < size; j++)
		{
			uw[i][j] /= buf;
			b[i][j] /= buf;
		}
	}
	for (i = 0; i < size; i++)
	{
		for (j = 0; j < size; j++)
		{
			uw[i][j] = b[i][j];
		}
	}
	return 1;
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
		case 'a': {d[i] = 't'; break; }
		case 't': {d[i] = 'a'; break; }
		case 'c': {d[i] = 'g'; break; }
		case 'g': {d[i] = 'c'; break; }
		default: d[i] = 'n';
		}
	}
	return 1;
}
int EvalMahControl(town* a, int nseq, int nseqb, int n_train, int n_cntrl, int* xporti, int* xportj, int*** seq, int*** seq_back, int olen, int* len, int* lenb, double** dav, double** dcv, char* filemat, FILE* outlog,FILE *out_roc, FILE *out_prc,char *file_for,double fp2, int iter, int irank)//char* filemap, 
{
	int i, k, n, m, o, b, u;
	double av[POPSIZE], buf[POPSIZE];
	double df[POPSIZE];

	for (k = 0; k < a->size; k++)
	{
		for (n = 0; n < a->size; n++)uw[k][n] = 0;
	}
	for (k = 0; k < a->size; k++)df[k] = 0;
	for (b = 0; b < n_train; b++)
	{
		m = xporti[b];		
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
			df[k] += fs[k];
		}
	}
	for (k = 0; k < a->size; k++)for (n = 0; n < a->size; n++)uw[k][n] /= n_train;
	for (k = 0; k < a->size; k++)df[k] /= n_train;
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
		av[k] = (df[k] + dav[a->tot[k].end - a->tot[k].sta][a->tot[k].num]) / 2;
		df[k] = df[k] - dav[a->tot[k].end - a->tot[k].sta][a->tot[k].num];
	}
	double* fp_rate;
	fp_rate = new double[n_cntrl];
	double* fp_count;
	fp_count = new double[n_cntrl];
	double* tp_sco;
	tp_sco = new double[n_cntrl];
	for (k = 0; k < n_cntrl; k++)fp_rate[k] = fp_count[k] = tp_sco[k] = 0;
	if (BackMat(a->size, outlog) == -1)
	{
		a->fit = 0;
		for (k = 0; k < n_cntrl; k++)
		{
			fp_count[k] = fp_rate[k] = 0.5;//n_cntrl_tot 
			tp_sco[k] = 0;
		}
		return 0;
	}
	double amahc = 0;
	for (k = 0; k < a->size; k++)
	{
		buf[k] = 0;
		for (n = 0; n < a->size; n++)buf[k] += uw[k][n] * df[n];
		amahc += buf[k] * df[k];
	}	
	for (k = 0; k < a->size; k++)
	{
		buf[k] /= amahc;
	}
	double sga_min = 0, sga_max = 0;
	for (k = 0; k < a->size; k++)
	{
		if (buf[k] < 0)sga_min += buf[k];
		else sga_max += buf[k];
	}
	double sga_raz = sga_max - sga_min;
	a->fprint_allfi_mat(filemat, olen, sga_min, sga_raz, buf);
	for (b = 0; b < n_cntrl; b++)fp_rate[b] = 0; 
	for (b = 0; b < n_cntrl; b++)fp_count[b] = 0;
/*	FILE* outmap;
	if ((outmap = fopen(filemap, "wt")) == NULL)
	{
		printf("Ouput file can't be opened!\n");
		exit(1);
	}
*/	int olen1 = olen - 1;
	for (b = 0; b < n_cntrl; b++)//control
	{
		u = xportj[b];
		int lenp = len[u] - olen1;
		double sco_pos = -1000;
//		int obest = 0, pbest = 0;
		for (m = 0; m < lenp; m++)
		{
			for (o = 0; o < 2; o++)
			{
				int gom = 1;
				double sco = 0;
				for (k = 0; k < a->size; k++)
				{
					double fs = 0;
					int rlenk = (a->tot[k].end - a->tot[k].sta + 1);
					for (n = a->tot[k].sta; n <= a->tot[k].end; n++)
					{
						int sym = seq[o][u][n + m];
						if (sym == -1)
						{
							gom = -1;
							break;
						}
						if (a->tot[k].num == sym)fs++;
					}
					if (gom == -1)break;
					if (fs > 0)
					{
						fs /= rlenk;
						sco += buf[k] * fs;
					}
				}
				if (gom == -1)sco = -1000;
				else sco = (sco - sga_min) / sga_raz;				
				if (sco > sco_pos)
				{
					sco_pos = sco;
	//				obest = o;
		//			pbest = m;
				}				
			}
		}		
		/*{
			char dir[] = "+-";
			char d[50];
			char letter[5] = "acgt";
			char err = 'n';
			memset(d, '\0', olen+1);			
			for (m = 0; m < olen1; m++)
			{
				int cod2 = seq[obest][u][pbest+m];
				int cod1 = cod2 / 4;
				if (cod1 >= 0 && cod1 <= 3)d[m] = letter[cod1];
				else d[m] = err;
			}
			{
				int cod2 = seq[obest][u][pbest + olen1 -1];
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
				pbest_sta = len[u] - pbest - olen1;
				pbest_end = len[u] - pbest;
			}
			fprintf(outmap, "%d\t%d\t%d\t%c\t%f\t%s\n", u + 1, pbest_sta, pbest_end, dir[obest], sco_pos, d);			
		}*/
		tp_sco[b] = sco_pos;
	}	
	//fclose(outmap);
	qsort(tp_sco, n_cntrl, sizeof(double), compare_qq2);
	int n_cntrl1 = n_cntrl - 1;
	int nseqn = 0;	
	for (b = 0; b < nseqb; b++)
	{
		int lenp = lenb[b] - olen;		
		double sco_max = -1000;
		for (m = 0; m <= lenp; m++)
		{
			double sco[2] = { 0,0 };
			int gom = 0;
			int mc;
			for (o = 0; o < 2; o++)
			{		
				if (o == 0)mc = m;
				else mc = lenp - m;// = 0 if m = max = lenp
				for (k = 0; k < a->size; k++)
				{
					double fs = 0;
					int rlenk = (a->tot[k].end - a->tot[k].sta + 1);
					for (n = a->tot[k].sta; n <= a->tot[k].end; n++)
					{
						int sym = seq_back[o][b][n + mc];
						if (sym == -1)
						{
							gom = -1;
							break;
						}
						if (a->tot[k].num == sym)fs++;
					}
					if (gom == -1)break;
					if (fs > 0)
					{
						fs /= rlenk;
						sco[o] += buf[k] * fs;
					}
				}
				if (gom == -1)break;
				sco[o] = (sco[o] - sga_min) / sga_raz;				
			}
			if (gom == 0)
			{
				nseqn++;
				double sco2 = Max(sco[0], sco[1]);
				if (sco_max < sco2)sco_max = sco2;
				if (sco2 >= tp_sco[n_cntrl1])
				{					
					for (k = n_cntrl1; k >=0; k--)
					{
						if (sco2 >= tp_sco[k])
						{
							fp_rate[k]++;							
						}
						else break;
					}
				}
			}

		}
		if (sco_max >= tp_sco[n_cntrl1])
		{			
			for (k = n_cntrl1; k >= 0; k--)
			{
				if (sco_max >= tp_sco[k])
				{
					fp_count[k]++;					
				}
				else break;
			}
		}
	}
	for (k = 0; k < n_cntrl; k++)
	{
		if (fp_rate[k] > 0)fp_rate[k] /= nseqn;
		else fp_rate[k] = 0.5 / (double)nseqn;
		fp_count[k] /= nseqb;		
	}
	fprintf(outlog, "Internal Size %d\n", a->size);
	for (k = 0; k < n_cntrl; k+=100)
	{
		fprintf(outlog,"ERR %f Count %d FPsites %f FPpeak %f TPpeak %f\n", tp_sco[k], k+1, fp_rate[k], fp_count[k],(double)(k+1)/n_cntrl);
	}
//	qsort(fp_rate, n_cntrl, sizeof(double), compare_qq);
//	qsort(fp_count, n_cntrl, sizeof(double), compare_qq);
//	qsort(tp_sco, n_cntrl, sizeof(double), compare_qq2);	
	double auc_roc=0, auc_pr =0;
	{
		fprintf(outlog, "\nROC\n");
		fprintf(out_roc, "\tROC_%d_%d_%d_%d\n", olen, a->size,iter+1,irank+1);
		fprintf(out_roc, "0\t0\n");
		double tproc_pred = 0;
		double fproc_pred = 0;
		for (n = 0; n < n_cntrl; n++)
		{
			if (fp_rate[n] > fproc_pred && (n == n_cntrl1 || fp_rate[n + 1] > fp_rate[n]))
			{
				double tproc_cur = (double)(n + 1) / n_cntrl;
				double fproc_cur = fp_rate[n];
				if (fproc_cur >= fp2)
				{
					double wei = (fp2 - fproc_pred) / (fproc_cur - fproc_pred);
					tproc_cur = tproc_pred + (tproc_cur - tproc_pred) * wei;
					fproc_cur = fp2;
				}
				double dauc = (tproc_cur + tproc_pred) * (fproc_cur - fproc_pred) / 2;
				fprintf(out_roc, "%g\t%f\n", fproc_cur, tproc_cur);
				auc_roc += dauc;
				if (fp_rate[n] >= fp2)break;
				tproc_pred = tproc_cur;
				fproc_pred = fproc_cur;
			}
		}
		auc_roc /= fp2;
		a->aucroc = auc_roc;
	}
	{		
		fprintf(out_prc, "\tPRC_%d_%d_%d_%d\n", olen, a->size, iter+1, irank+1);
		double tpc, tpc_pred = 0;
		double fpc, fpc_pred = 0;
		double prec_pred = 1;
		int prc_count = 1;
		double fp_pred = 0;
		double prec_exp = 0.5;
		fprintf(out_prc, "0\t1\n");
		int n_cntrl1 = n_cntrl - 1;
		for (i = 0; i < n_cntrl; i++)
		{
			int i1 = i + 1;
			if (i == n_cntrl1 || (fp_rate[i1] != fp_rate[i] || tp_sco[i1] != tp_sco[i]))
			{
				tpc = (double)(i1) / n_cntrl;
				fpc = fp_count[i];
				if (fp_rate[i] >= fp2)
				{
					double dtpc = tpc - tpc_pred;
					double dfpc = fpc - fpc_pred;
					double wei = (fp2 - fp_pred) / (fp_rate[i] - fp_pred);
					tpc = tpc_pred + wei * dtpc;
					fpc = fpc_pred + wei * dfpc;
				}
				double prec_cur = tpc / (tpc + fpc);
				double prec_av = (prec_pred + prec_cur) / 2;
				double dauc = (prec_av - prec_exp) * (tpc - tpc_pred);
				auc_pr += dauc;
				fprintf(out_prc, "%f\t%f\n", tpc, prec_cur);
				tpc_pred = tpc;
				fpc_pred = fpc;
				fp_pred = fp_rate[i];
				prec_pred = prec_cur;
				prc_count++;
				if (fp_rate[i] >= fp2)break;
			}
		}
		auc_pr *= 2;
		a->auprc = auc_pr;
	}
	delete[] fp_rate;
	delete[] fp_count;
	delete[] tp_sco;
	return 1;
}
double EvalMahFIT(town* a, int n_train, int octa, int infc, int* xporti, int*** seq, int olen, double** dav, double** dcv, double** octa_prow, double *pexp, int* len, FILE* outlog)//double *hoxa_wei, qbs *qps, double **frp, double *octa_rat,
{
	int k, n, m, b;
	double f1[POPSIZE], df[POPSIZE], f2[POPSIZE];
	double pfm[MOTLEN][16];

	int olen1 = olen - 1;
	for (k = 0; k < olen1; k++)for (n = 0; n < 16; n++)pfm[k][n] = 0;	

	//fprintf(out,"Enter Fit1 %d\n",n_train);
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
	//fprintf(out,"Enter Fit2 %d\n", n_train);
	for (b = 0; b < n_train; b++)
	{
		m = xporti[b];
		//if (len[m] < olen)continue;
		int ori = a->ori[m];
		int pos = a->pos[m];
		for (k = 0; k < olen1; k++)
		{
			int let = seq[ori][m][k + pos];
			if (let != -1)pfm[k][let]++;
		}
		double fs[POPSIZE];
		//fprintf(out,"Seq %d\t",b);
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
			//fprintf(out,"%d %g\t", k, frp[b][k]);
		}
		//fprintf(out,"\n");		
		for (k = 0; k < a->size; k++)
		{
			for (n = 0; n < a->size; n++)
			{
				uw[k][n] += fs[k] * fs[n];
			}
			f1[k] += fs[k];
		}
	}	
	//exit(1);
	for (k = 0; k < olen1; k++)
	{
		int  n_train1 = n_train + 1;
		for (n = 0; n < 16; n++)
		{
			pfm[k][n] += pexp[n];
			/*if (pfm[k][n] == 1)
			{
				int yy = 0;
			}*/
			pfm[k][n] /= n_train1;
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
				a->inf += pfm[k][n] * log(pfm[k][n] / pexp[n]);
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
	//fprintf(out,"Av1,2\n");
	for (k = 0; k < a->size; k++)
	{
		df[k] = f1[k] - f2[k];
		//fprintf(out,"%g %g\t", av[k], df[k]);
	}
	if (BackMat(a->size, outlog) == -1) { a->fit = 0; return 0; }
	a->mah = 0;
	for (k = 0; k < a->size; k++)
	{
		double buf = 0;
		for (n = 0; n < a->size; n++)buf += uw[k][n] * df[n];
		a->mah += buf * df[k];
	}
	double oct = 0;
	int sum = 0;
	for (b = 0; b < n_train; b++)
	{
		//m = qps[b].n;
		m = xporti[b];
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
		oct += wei;
		sum++;
		//fprintf(out,"oct %f sco %f peak %d cep %d pos %d pro %f wei %g\n", oct, qps[b].q, m, a->ori[m], a->pos[m], octa_prow[a->ori[m]][m][a->pos[m]], hoxa_wei[b]);
	}
	oct /= sum;
	oct = pow(10, oct);
	a->fpr = oct;
	a->fit = a->mah * a->fpr * a->inf;
	return a->fit;
}
int MutOlig0(town* a)
{
	int cy = 0, r1, i, gom;
	int num1;
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
	//	fprintf(out,"%d%d ",a->odg[w0],a->odg[w]);
	a->odg[w0]--;
	a->odg[w]++;
	//	fprintf(out,"%d%d \n",a->odg[w0],a->odg[w]);
	a->tot[n1].sta = a2;
	a->tot[n1].end = b2;
	return n1;
}
int MutRegShift(town* a, int n_train, int* len, int* xporti, int olen, int& npeak, int& nori, int& npos)
{
	//fprintf(out,"In Peak %d Ori %d Pos %d ",npeak,nori,npos);
	int r1, r2, cy = 0;
	r1 = rand() % n_train;
	r1 = xporti[r1];
	npeak = r1;
	int lenr = len[r1] - olen + 1;
	r2 = rand() % lenr;
	npos = r2;
	if (r2 == a->pos[r1])
	{
		nori = a->ori[r1] = 1 - a->ori[r1];
		//fprintf(out,"Out1 Peak %d Pos %d Ori %d",npeak, npos,nori);
		return 1;
	}
	else
	{
		a->pos[r1] = r2;
		int r3 = rand() % 2;
		if (r3 == 1)nori = a->ori[r1] = 1 - a->ori[r1];
		else nori = a->ori[r1];
		//fprintf(out,"Out2 Peak %d Pos %d Ori %d",npeak, npos,nori);
		return 1;
	}
}
/*int MutRegShiftHoxa(town *a, int n_train, int *xporti, int &npeak, int &nori, int &npos, int *len_octa, int **octa_prowb, int *len, int olen)
{
	//fprintf(out,"In Peak %d Ori %d Pos %d ",npeak,nori,npos);
	int r1, r2, oln1 = olen -1;
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
	//fprintf(out,"Out2 Peak %d Pos %d Ori %d",npeak, npos,nori);
	return 1;

}*/
int MutRegShiftHoxaW(town* a, int n_train, int* xporti, int& npeak, int& nori, int& npos, int* len_octa, int** octa_prowb, int** octa_prows)
{
	//fprintf(out,"In Peak %d Ori %d Pos %d ",npeak,nori,npos);
	int r1, r2, j;
	r1 = rand() % n_train;
	npeak = xporti[r1];
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
	//fprintf(out,"Out2 Peak %d Pos %d Ori %d",npeak, npos,nori);
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
		//		fprintf(out,"\n");
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
int Reco2Peak(town* a1, town* a2, int n_train, int* xporti)
{
	int r1, r2 = 0;
	r1 = rand() % n_train;
	do
	{
		r2 = rand() % n_train;
	} while (r1 == r2);
	r1 = xporti[r1];
	r2 = xporti[r2];
	if (a1->pos[r1] == a2->pos[r2] && a1->ori[r1] == a2->ori[r2])return -1;
	MixI(&a1->ori[r1], &a2->ori[r1]);
	MixI(&a1->pos[r1], &a2->pos[r1]);
	//	MixI(&a1->posc[r1], &a2->posc[r1]);
	MixI(&a1->ori[r2], &a2->ori[r2]);
	MixI(&a1->pos[r2], &a2->pos[r2]);
	//	MixI(&a1->posc[r2], &a2->posc[r2]);
	return 1;
}
void MixPop(town* a, town* b)
{
	town c = *a;
	*a = *b;
	*b = c;
}
int GomTown(town a, town b, int nseq, int* xporti, int check_lpd)
{
	int i, j;
	for (i = 0; i < nseq; i++)
	{
		j = xporti[i];
		if (b.pos[j] != a.pos[j])return 0;
		if (b.ori[j] != a.ori[j])return 0;
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
	memset(d, 0, sizeof(d));
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
	memset(d, 0, sizeof(d));
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
	for (i = 0; i < 2; i++)memset(d[i], 0, sizeof(d[i]));
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
				/*for (j = 0; j < 2; j++)
				{
					strcpy(peak_real[j][n], d[j]);
					peak_real[j][n][len[n]] = '\0';
				}*/
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
							//	fprintf(out,"ReadSeq! Input file %s peak %d strand %d position %d error!\n%s\n", file, n, j, i, d[j]);
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
				if (check == -1)
				{
					printf("ReadSeq! Unusual symbol, peak %d ignored\n%s\n", nn + 1, d[0]);
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
void ReadSeqBack(char* file, int nseq, int* len, int*** seq_back, int olen, int reg_max, double** dav, double** dcv, double* octaf, double *pexp, int* octa1, int octa, int octa_size, int len_peak_max, FILE* outlog)
{
	char l[SEQLEN], d[SEQLEN], head[400];
	int m, i, j, k, v, fl = 0;
	int n_backr[LPDLEN];
	double freq[16];
	FILE* in;

	if ((in = fopen(file, "rt")) == NULL)
	{
		printf("ReadSeq! Input file %s can't be opened!\n", file);
		exit(1);
	}
	char symbol = fgetc(in);
	rewind(in);
	for (j = 0; j < reg_max; j++)n_backr[j] = 0;
	int nn = 0, n = 0;
	memset(d, 0, sizeof(d));
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
				for (j = 0; j < reg_max; j++)n_backr[j] += (len[n] - j + 1);
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
						if (cod[0] >= 0 && cod[1] >= 0)code = 4 * cod[0] + cod[1];
						//	else
						//	{
							//	fprintf(out,"ReadSeqBack! Input file %s peak %d strand %d position %d error!\n%s\n", file, n, m, k, d);
								//exit(1);
							//}
						seq_back[m][n][k] = code;
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
								int sb = seq_back[m][n][k + v];
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
					printf("ReadSeqBack!Unusual symbol, peak %d ignored\n%s\n", nn + 1, d);
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
int TestArgv(char* arg1, char* arg2, char* arg3, char* arg4, char* arg5, char* arg6, char* arg7, char* arg8, char* arg9, char* arg10, char* arg11, char* arg12,
	char* val1, char* val2, char* val3, int val4, int val5, int val6, int val7, double val8, int val9, char* val10, int val11, char* val12, FILE *out)
{
	int ret = 0;
	if (strcmp(arg1, val1) != 0)
	{
		fprintf(out, "%s is not %s\n", arg1, val1);
		return 1;
	}
	if (strcmp(arg2, val2) != 0)
	{
		fprintf(out, "%s is not %s\n", arg2, val2);
		return 2;
	}
	if (strcmp(arg3, val3) != 0)
	{
		fprintf(out, "%s is not %s\n", arg3, val3);
		return 3;
	}
	int v4 = atoi(arg4);
	if (v4 != val4)
	{
		fprintf(out, "%d is not %d\n", v4, val4);
		return 4;
	}
	int v5 = atoi(arg5);
	if (v5!=val5)
	{
		fprintf(out, "%d is not %d\n", v5, val5);
		return 5;
	}
	int v6 = atoi(arg6);
	if (v6 != val6)
	{
		fprintf(out, "%d is not %d\n", v6, val6);
		return 6;
	}
	int v7 = atoi(arg7);
	if (v7 != val7)
	{
		fprintf(out, "%d is not %d\n", v7, val7);
		return 7;
	}
	double v8 = atof(arg8);
	if (v8 != val8)
	{
		fprintf(out, "%f is not %f\n", v8, val8);
		return 8;
	}
	int v9 = atoi(arg9);
	if (v9 != val9)
	{
		fprintf(out, "%d is not %d\n", v9, val9);
		return 9;
	}
	if (strcmp(arg10, val10) != 0)
	{
		fprintf(out, "%s is not %s\n", arg10, val10);
		return 10;
	}
	int v11 = atoi(arg11);
	if (v11 != val11)
	{
		fprintf(out, "%d is not %d\n", v11, val11);
		return 11;
	}
	if (strcmp(arg12, val12) != 0)
	{
		fprintf(out, "%s is not %s\n", arg12, val12);
		return 13;
	}
	return 0;
}
int main(int argc, char* argv[])
{
	int* len, nseq, nseqb, * lenb, i, j, k, n, m;
	char file_for[500], file_back[500], path_fasta[500], path_out[500], pfile_for[500], pfile_back[500], pfile_log[500], file_log[500];
	char file_out_roc[500], file_out_prc[500], file_out_auc[500];
	char filemat[500];// , filemap[500];
	int*** seq_real, *** seq_back;
	double** dav;//dinucl.content background
	double** dcv;//self covariations for regions LPD
	int** octa_prowb, * len_octa, ** octa_prows;// octa position lists, octa position counts, weight sums
	double** octa_pro1, ** octa_prow, * thr_octa;// , *hoxa_wei;
	double pexp[16];

	FILE* outlog, * out_roc, * out_prc, * out_auc;

	if (argc != 14)
	{
		puts("Sintax: 1path_both_fasta 2,3files_forground,background 4int LPD_length 5,6int no. of LPD_min,dif 7int motif_len 8double FPR threshold 9int olig_background 10int infc(0no,1yes) 11path_out 12int max_peak_len 13file log");//  5<pop_size>
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
	int size_start = atoi(argv[5]);// no. of LPD min	
	int size_dif = atoi(argv[6]);// no. of LPD dif
	int size_end = size_start + (CELL - 1) * size_dif;
	int olen = atoi(argv[7]);// dlina motiva
	double fp2 = atof(argv[8]);// 0.01 = 1 site per 100 base pairs, FPR threshold for pAUC	
	double ratio_train_to_control = -1; // atof(argv[8]); // value -1 means equal sizes of training and control subsets, odd/even peaks are used either for training and control subsets
	int iteration = 2; // atoi(argv[9]);//total no. of jack-knife test			
	int octa = atoi(argv[9]);
	int infc = atoi(argv[10]); // information content in GA fitness
	strcpy(path_out, argv[11]);
	int size0[CELL];
	for (i = 0; i < CELL; i++)size0[i] = size_start + i * size_dif;	
	int len_peak_max = atoi(argv[12]); //2500;
	strcpy(pfile_log, path_out);
	strcpy(file_log, argv[13]);
	strcat(pfile_log, file_log);
	if ((outlog = fopen(pfile_log, "wt")) == NULL)
	{
		fprintf(outlog, "Input file %s can't be opened!\n", pfile_log);
		exit(1);
	}
	time_t tnow = time(NULL);
	fprintf(outlog, "%s", ctime(&tnow));
	fprintf(outlog, "Start1\n");
	{
		int testar = TestArgv(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], argv[8], argv[9], argv[10], argv[11], argv[12],
			path_fasta, file_for, file_back, reg_max, size_start, size_dif, olen, fp2, octa, path_out, len_peak_max, file_log, outlog);
		if (testar != 0)
		{
			printf("TestArgv Start point: %s argument %d error\n", argv[testar],testar);
			exit(1);
		}
	}
	srand((unsigned)time(NULL));
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
	nseq = nseqb = 0;
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
	//
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
	seq_back = new int** [2];
	if (seq_back == NULL) { puts("Out of memory..."); exit(1); }
	for (i = 0; i < 2; i++)
	{
		seq_back[i] = new int* [nseqb];
		for (j = 0; j < nseqb; j++)
		{
			seq_back[i][j] = new int[lenb[j] - 1];
			if (seq_back[i][j] == NULL) { puts("Out of memory..."); exit(1); }
		}
	}
	for (i = 0; i < 2; i++)
	{
		for (j = 0; j < nseq; j++)for (k = 0; k < len[j] - 1; k++)seq_real[i][j][k] = -1;
		for (j = 0; j < nseqb; j++)for (k = 0; k < lenb[j] - 1; k++)seq_back[i][j][k] = -1;
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
	ReadSeq(pfile_for, nseq, len, seq_real, olen, octa_av, octa1, octa, octa_size, octa_pro1, len_peak_max, outlog);//peak_real, 
	for (i = 0; i < octa_size; i++)octa_rat[i] = log10(octa_av[i]);
	for (i = 0; i < octa_size; i++)octa_av[i] = 0;
	for (i = 0; i < octa_size; i++)octa1[i] = 0;
	for (i = 0; i < 16; i++)pexp[i] = 0;
	ReadSeqBack(pfile_back, nseqb, lenb, seq_back, olen, reg_max, dav, dcv, octa_av, pexp, octa1, octa, octa_size, len_peak_max, outlog);
	for (i = 0; i < octa_size; i++)octa_rat[i] -= log10(octa_av[i]);
	for (i = 0; i < nseq; i++)
	{
		int leni = len[i] - octa + 1;
		for (n = 0; n < leni; n++)
		{
			int i_sost = (int)octa_pro1[i][n];
			octa_pro1[i][n] = octa_rat[i_sost];
		}
	}
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
	if (ratio_train_to_control == 0)
	{
		iteration = Min(nseq, iteration);
	}
	int* n_train;
	n_train = new int[iteration];
	if (n_train == NULL) { puts("Out of memory..."); exit(1); }
	int* n_cntrl;
	n_cntrl = new int[iteration];
	if (n_cntrl == NULL) { puts("Out of memory..."); exit(1); }
	int iter;
	if (ratio_train_to_control == 0) //Leave-one-out cross-validation
	{
		int nseq1 = nseq - 1;
		for (i = 0; i < iteration; i++) { n_train[i] = nseq1; n_cntrl[i] = 1; }
	}
	else
	{
		if (ratio_train_to_control > 0)// random splits cross-validation
		{
			int n_cnt = (int)(nseq / ratio_train_to_control);
			int n_trn = nseq - n_cnt;
			for (iter = 0; iter < iteration; iter++)
			{
				n_train[iter] = n_trn;
				n_cntrl[iter] = n_cnt;
			}
		}
		else //two-fold cross-validation
		{
			int n_less = nseq / 2;
			int n_more = nseq - n_less;
			for (iter = 0; iter < iteration; iter++)
			{
				if (iter % 2 == 0)
				{
					n_train[iter] = n_less;
					n_cntrl[iter] = n_more;
				}
				else
				{
					n_train[iter] = n_more;
					n_cntrl[iter] = n_less;
				}
			}
		}
	}
	int n_cnt_tot = 0;
	for (i = 0; i < iteration; i++)n_cnt_tot += n_cntrl[i];
	int n_trn_tot = 0;
	for (i = 0; i < iteration; i++)n_trn_tot += n_train[i];
	int n_both_sam = 0;
	for (i = 0; i < iteration; i++)
	{
		n_both_sam += n_cntrl[i];
		n_both_sam += nseqb;
	}

	int* xport;
	xport = new int[nseq];
	if (xport == NULL) { puts("Out of memory..."); exit(1); }
	int* xporti;//train
	xporti = new int[nseq];
	if (xporti == NULL) { puts("Out of memory..."); exit(1); }
	int* xportj;//control
	xportj = new int[nseq];
	if (xportj == NULL) { puts("Out of memory..."); exit(1); }
	{
		char word[] = "acgt";
		GetWords(2, 0, 16, word);
	}	
	pop = new town * *[iteration];
	if (pop == NULL) { puts("Out of memory...!"); exit(1); }
	for (iter = 0; iter < iteration; iter++)
	{
		pop[iter] = new town * [CELL];
		if (pop[iter] == NULL) { puts("Out of memory...!"); exit(1); }
		for (i = 0; i < CELL; i++)
		{
			pop[iter][i] = new town[MEGE];
			if (pop[iter][i] == NULL) { puts("Out of memory...!"); exit(1); }
		}
	}
	for (iter = 0; iter < iteration; iter++)
	{
		for (i = 0; i < CELL; i++)
		{
			for (j = 0; j < MEGE; j++)
			{
				pop[iter][i][j].mem_in(nseq);
			}
		}
	}
	for (i = 0; i < 2; i++)
	{
		det2[i].mem_in(nseq);
	}
	det1.mem_in(nseq);
	{
		char name[200];
		memset(name, '\0', sizeof(name));
		for (i = 0;; i++)
		{
			if (file_for[i] == '.') { name[i] = '\0'; break; }
			if (file_for[i] == '\n') { name[i] = '\0'; break; }
			if (file_for[i] == '\0') { name[i] = '\0'; break; }
			name[i] = file_for[i];
		}
		char add_auc[500], add_roc1[500], add_prc1[500]; //add_prc[500],
		strcpy(add_auc, "_auc_bs.txt"); // common
		strcpy(add_prc1, "_prc_bs1.txt"); //PRC
		strcpy(add_roc1, "_roc_bs1.txt"); //ROC
		memset(file_out_roc, 0, sizeof(file_out_roc));
		strcpy(file_out_roc, path_out);
		//strcat(file_out_roc, argv[2]);
		strcat(file_out_roc, name);
		strcat(file_out_roc, "_");
		strcat(file_out_roc, argv[7]);
		strcat(file_out_roc, add_roc1);
		memset(file_out_prc, 0, sizeof(file_out_prc));
		strcpy(file_out_prc, path_out);
		//strcat(file_out_prc, argv[2]);
		strcat(file_out_prc, name);
		strcat(file_out_prc, "_");
		strcat(file_out_prc, argv[7]);
		strcat(file_out_prc, add_prc1);
		memset(file_out_auc, 0, sizeof(file_out_auc));
		strcpy(file_out_auc, path_out);
		//strcat(file_out_auc, argv[2]);
		strcat(file_out_auc, name);
		strcat(file_out_auc, "_");
		strcat(file_out_auc, argv[7]);
		strcat(file_out_auc, add_auc);		
	}
	if ((out_roc = fopen(file_out_roc, "wt")) == NULL)
	{
		printf("Output file can't be opened!\n");
		exit(1);
	}
	if ((out_prc = fopen(file_out_prc, "wt")) == NULL)
	{
		printf("Output file can't be opened!\n");
		exit(1);
	}
	if ((out_auc = fopen(file_out_auc, "wt")) == NULL)
	{
		printf("Output file can't be opened!\n");
		exit(1);
	}
	int n_train_max = 0;
	for (iter = 0; iter < iteration; iter++)
	{
		if (n_train[iter] > n_train_max)n_train_max = n_train[iter];
	}
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
	for (i = 0; i < nseq; i++)for (n = 0; n < len[i]; n++)octa_prowb[i][n] = -1;
	octa_prows = new int* [nseq];
	if (octa_prows == NULL) return -1;
	for (i = 0; i < nseq; i++)
	{
		octa_prows[i] = new int[len[i]];
		//octa_prows[i] = new int[len[i]];
		if (octa_prows[i] == NULL) return -1;
	}
	int pair_all_max = 0;
	int** pair_d; // pair 1 & 2 population ranks of individuals in each cell				
	int* pair_take;
	int* first_cell, * second_cell;
	double jwei;
	double jwei0[3];// = { 1.1, 1.2, 1.3 };
	{
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
		fprintf(outlog, "RecWei\t");
		for (j = jmax; j >= 0; j--)fprintf(outlog, "%d %d\t", j, rwei[j]);
		fprintf(outlog, "\n");
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
	first_cell = new int[pair_all_max];
	if (first_cell == NULL) return -1;
	second_cell = new int[pair_all_max];
	if (second_cell == NULL) return -1;
	pair_take = new int[pair_all_max];
	if (pair_take == NULL) return -1;
	pair_d = new int* [pair_all_max];
	if (pair_d == NULL) return -1;
	for (m = 0; m < pair_all_max; m++)
	{
		pair_d[m] = new int[2];
		if (pair_d[m] == NULL) return -1;
	}
	{
		int testar = TestArgv(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], argv[8], argv[9], argv[10], argv[11], argv[12],
			path_fasta, file_for, file_back, reg_max, size_start, size_dif, olen, fp2, octa, path_out, len_peak_max, file_log, outlog);
		if (testar != 0)
		{
			printf("TestArgv Enter Iterations: %s argument %d error\n", argv[testar], testar);
			exit(1);
		}
	}
	{
		int len_tot = 0, len_wei = 0;
		for (i = 0; i < nseq; i++)
		{
			int leni = len[i] - olen + 1;
			len_tot += leni;
			int odif = olen - octa + 1;
			for (n = 0; n < leni; n++)
			{
				octa_prow[i][n] = 0;
				for (k = 0; k < odif; k++)octa_prow[i][n] += octa_pro1[i][n + k];
				octa_prow[i][n] /= odif;
			}
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
					//sumw += dw;					
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
		fprintf(outlog, "FractHoxa %f\n", (double)len_wei / len_tot);
		//double auc_roc_len = 0, auc_pr_len = 0;
		{
			int big_exit1 = 1;// local exit (separ +-) global exit (separation do not exceeded the previous run)
			double fit_prev[CELL], fit_after_mut[CELL];
			for (iter = 0; iter < iteration; iter++)
			{
				fprintf(outlog, "\n%s vs %s\tWindowlength %d\tBootStrap %d  the last %d\t", file_for, file_back, olen, iter + 1, iteration);
				fprintf(outlog, "N LPD [%d", size_start);
				for (i = 1; i < CELL; i++)fprintf(outlog, ",%d", size0[i]);
				fprintf(outlog, "]\tDeg %d\tEli %d\tCELL %d\n", MEGE, MEGE, CELL);
				if (ratio_train_to_control > 0)// random splits cross-validation
				{
					for (k = 0; k < n_cntrl[iter]; k++)xport[k] = 0;
					for (k = n_cntrl[iter]; k < nseq; k++)xport[k] = 1;
					BigMixI(xport, nseq);
				}
				else
				{
					if (ratio_train_to_control == 0) //Leave-one-out cross-validation
					{
						for (k = 0; k < nseq; k++)xport[k] = 1;//train
						xport[iter] = 0;//control
					}
					else //two-fold cross-validation
					{

						if (iter % 2 == 0)
						{
							xport[0] = 0;//control
						}
						else
						{
							xport[0] = 1;//train
						}
						for (k = 1; k < nseq; k++)
						{
							xport[k] = 1 - xport[k - 1];
						}
					}
				}
				{
					int it = 0;//train
					for (k = 0; k < nseq; k++)
					{
						if (xport[k] == 1)xporti[it++] = k;
					}
				}
				//initiation		
				int gen = 0;
				int success_m;
				int success_r1[NREC];				
				int check_lpd = 1;
				int stop_oi[CELL][MEGE], stop_li[CELL][MEGE], stop_pi[CELL][MEGE];
				int success_rcell[CELL];
				for (i = 0; i < CELL; i++)success_rcell[i] = 0;
				int success_ri[CELL][MEGE];
				for (i = 0; i < CELL; i++)for (j = 0; j < MEGE; j++)stop_oi[i][j] = stop_li[i][j] = stop_pi[i][j] = 0;
				for (i = 0; i < CELL; i++)for (j = 0; j < MEGE; j++)success_ri[i][j] = 1;
				//int restart_full = 0, restart_half = 0;
				do
				{
					int success_o, success_l, success_p;
					if (big_exit1 == 1)
					{
						for (i = 0; i < CELL; i++)
						{
							for (j = 0; j < MEGE; j++)
							{
								int err = 1;
								do
								{
									err = 1;
									int gom = 0;
									det1.init_rand_hoxa(nseq, olen, size0[i], reg_max, len_octa, len, octa_prowb, octa_prows);
									if (det1.check(0, reg_max, outlog) == -1)
									{
										det1.check(0, reg_max, outlog);
										fprintf(outlog, "Population error!\n");
										det1.print_all(reg_max, nseq, outlog);
										exit(1);
									}
									else
									{
										for (m = 0; m < j; m++)
										{
											gom = GomTown(det1, pop[iter][i][m], n_train[iter], xporti, 1);
											if (gom == -1)
											{
												break;
											}
										}
									}
									if (gom != -1)
									{
										EvalMahFIT(&det1, n_train[iter], octa, infc, xporti, seq_real, olen, dav, dcv, octa_prow, pexp, len, outlog);//hoxa_wei,qps, frp, octa_rat,
										det1.get_copy(&pop[iter][i][j], nseq, reg_max);
										err = 0;
									}
								} while (err == 1);
							}
						}
						big_exit1 = 0;
						for (i = 0; i < CELL; i++)qsort((void*)(&pop[iter][i][0]), MEGE, sizeof(pop[iter][i][0]), compare_pop);
						for (i = 0; i < CELL; i++)pop[iter][i][0].print_sta(reg_max, outlog);
					}
					for (i = 0; i < CELL; i++)fit_prev[i] = pop[iter][i][0].fit;
					success_o = success_l = success_p = success_m = 0;
					double ratio2_gen0 = 0.01, ratio_rec_cycle, ratio_mut_cycle = 0.0001;
					int step, step_max, step_max_tot = 0;
					int mege_rec;
					int kn_train_max = Max(n_train_max, 500);
					int sr_step = 400 * kn_train_max;//200000 if n_train_max = 500		;
					int n_rec_cycle_max = 4000 * kn_train_max;//2000000 if n_train_max = 500							
					if (gen == 0)
					{
						step = 2000;
						step_max = 20 * kn_train_max; //10000 if n_train_max = 500						
						mege_rec = MEGE / 5;
						jwei = jwei0[0];
						ratio_rec_cycle = 0.0025;
						sr_step = 400 * kn_train_max;
					}
					else
					{
						if (gen == 1)
						{
							jwei = jwei0[1];
							step = 4000;//25000;
							step_max = 200 * kn_train_max;//200 -> 100000, 500 -> 250000 if n_train_max = 500							
							mege_rec = MEGE / 5;
							ratio_rec_cycle = 0.001;
							sr_step = 1000 * kn_train_max;
							//for (i = 0; i < CELL; i++)for (j = 0; j < MEGE; j++)stop_oi[i][j] = stop_li[i][j] = stop_pi[i][j] = 0;
						}
						else
						{
							jwei = jwei0[2];
							step = 5000;
							step_max = 200 * kn_train_max;//200 -> 100000, 1000 -> 500000 if n_train_max = 500							
							mege_rec = MEGE / 5;														
							ratio_rec_cycle = 0.001;
							sr_step = 1000 * kn_train_max;
							//for (i = 0; i < CELL; i++)for (j = 0; j < MEGE; j++)stop_li[i][j] = stop_pi[i][j] = 0;//stop_oi[i][j] = 
						}
					}
					if (step_max < step)step_max = step;
					double ratio_thr = 1 / (double)step;
					double ratio_thr_r[2];
					for (i = 0; i < 2; i++)ratio_thr_r[i] = ratio_thr;
					int step2 = 2 * step;					
					int n_mut_tot = 0, success_m_tot = 0, mdo = 1;
					double sm0_rate = 0.001;
					int mut_jump[CELL];
					for (j = 0; j < CELL; j++)mut_jump[j] = 0;
					double exp_rec_rate[CELL][MEGE];
					for (j = 0; j < CELL; j++)for (i = 0; i < MEGE; i++)exp_rec_rate[j][i] = 0;
					double fit_pre_rec[CELL][MEGE];
					//mutations																													
					//if(gen<=1)
					{
						fprintf(outlog, "Enter mutations Gen %d Iteration %d\n", gen + 1, iter + 1);
						{
							int testar = TestArgv(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], argv[8], argv[9], argv[10], argv[11], argv[12],
								path_fasta, file_for, file_back, reg_max, size_start, size_dif, olen, fp2, octa, path_out, len_peak_max, file_log, outlog);
							if (testar != 0)
							{
								printf("TestArgv Before Mut: %s argument %d error\n", argv[testar], testar);
								exit(1);
							}
						}
						int asuccess[NMUT], atry[NMUT];
						for (k = 0; k < NMUT; k++)asuccess[k] = atry[k] = 0;
						int success_mi[CELL][MEGE], try_mi[CELL][MEGE];
						for (j = 0; j < CELL; j++)for (i = 0; i < MEGE; i++)success_mi[j][i] = try_mi[j][i] = 0;
						int mdo_cell[CELL];
						for (j = 0; j < CELL; j++)mdo_cell[j] = 1;
						do
						{
							step_max_tot += step_max;
							for (j = 0; j < CELL; j++)for (i = 0; i < MEGE; i++)pop[iter][j][i].sm = pop[iter][j][i].tm = 0;
							if (gen == 0)fprintf(outlog, "Mut cycle %d\n", step_max_tot / step_max);
							for (j = 0; j < CELL; j++)
							{
								if (mdo_cell[j] != 1)continue;
								for (i = 0; i < MEGE; i++)
								{
									if (stop_pi[j][i] == 0 || (stop_li[j][i] == 0 || stop_oi[j][i] == 0))
									{										
										int success_o_local = 0, success_l_local = 0, success_p_local = 0;
										int success_m_local = 0;
										double ratio[NMUT];
										if (stop_oi[j][i] == 0)ratio[0] = 1;
										else ratio[0] = 0;
										if (stop_li[j][i] == 0)ratio[1] = 1;
										else ratio[1] = 0;
										if (stop_pi[j][i] == 0)ratio[2] = 1;
										else ratio[2] = 0;
										int step_success[NMUT];
										int step_try[NMUT];
										for (k = 0; k < NMUT; k++)step_success[k] = step_try[k] = 0;
										int n_mut_here = 0;
										double fit_mut_prev0 = pop[iter][j][i].fit;
										int n_mut_min_cyc = step;
										int m_iter = 0;
										while (ratio[2] > ratio_thr || (ratio[0] > ratio_thr || ratio[1] > ratio_thr))
										{
											pop[iter][j][i].get_copy(&det1, nseq, reg_max);
											int sm;											
											if (stop_oi[j][i] == 1)
											{
												if (stop_li[j][i] == 1)sm = 2;
												else
												{
													sm = rand() % 2;// 1 or 2
													sm++;
												}
											}
											else
											{												
												if (stop_li[j][i] == 1)
												{
													sm = rand() % 2;
													if (sm == 1)sm++;// 0 or 2
												}
												else sm = rand() % 3;
											}
											/*{
												int rs = step + step_success[2];
												if (ratio[1] > ratio_thr && stop_li[j][i] == 0)rs += step + step_success[1];
												if (ratio[0] > ratio_thr && stop_oi[j][i] == 0)rs += step + step_success[0];
												int rs1 = rand() % rs;
												if (rs1 < step + step_success[2])sm = 2;
												else
												{
													if (rs1 < step2 + step_success[2] + step_success[1])sm = 1;
													else sm = 0;
												}
											}*/
										//	det1.print_all2(reg_max);
											int muto;
											int npeak = -1, nori = -1, npos = -1;
											if (sm == 0)muto = MutOlig0(&det1);
											else
											{
												if (sm == 1)muto = MutCry0(&det1, olen, reg_max);
												else muto = MutRegShiftHoxaW(&det1, n_train[iter], xporti, npeak, nori, npos, len_octa, octa_prowb, octa_prows);
											}
											int gom = 0;
											if (muto != -1)
											{
												int muto1 = det1.order(muto);
												if (muto1 == -1)
												{
													puts("MUTO Error");
													gom = -1;
												}
												if (gom == 0)
												{
													for (m = 0; m < MEGE; m++)
													{
														if (m == i)continue;
														gom = GomTown(det1, pop[iter][j][m], n_train[iter], xporti, check_lpd);
														if (gom == -1) { break; }
													}
												}
											}
											double ret = 0, dd = 0;											
											if (gom != -1)
											{												
												ret = EvalMahFIT(&det1, n_train[iter], octa, infc, xporti, seq_real, olen, dav, dcv, octa_prow, pexp, len, outlog);//hoxa_wei,qps, frp, octa_rat,
												if (ret > 0)
												{
													n_mut_here++;
													step_try[sm]++;
													try_mi[j][i]++;
													dd = det1.fit / pop[iter][j][i].fit;
													pop[iter][j][i].tm++;
													if (dd > 1)
													{
														/*if (det1.fit > pop[iter][j][0].fit)
														{
															printf("Cell %d Rank %d SM %d DD %f\n", j + 1, i + 1, sm, dd);
															det1.print_all2(reg_max);
															pop[iter][j][i].print_all2(reg_max);
														}*/
														step_success[sm]++;
														pop[iter][j][i].sm++;
														if (sm == 0)success_o_local++;
														else
														{
															if (sm == 1)success_l_local++;
															else success_p_local++;
														}
														success_m_local++;
														det1.get_copy(&pop[iter][j][i], nseq, reg_max);
														success_mi[j][i]++;
													}
												}
												//else pop[iter][j][i].get_copy(&det1, nseq, reg_max);
											}
											//else pop[iter][j][i].get_copy(&det1, nseq, reg_max);											
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
												if (ratio[0] <= ratio_thr)stop_oi[j][i] = 1;
												if (ratio[1] <= ratio_thr)stop_li[j][i] = 1;
												if (ratio[2] <= ratio_thr)stop_pi[j][i] = 1;
												double ratio_per_cycle = pop[iter][j][i].fit / fit_mut_prev0 - 1;
												for (k = 0; k < 3; k++)atry[k] += step_try[k];
												for (k = 0; k < 3; k++)asuccess[k] += step_success[k];
												int step_success_tot = step_success[0] + step_success[1] + step_success[2];
												if (step_success_tot < n_mut_min_cyc)n_mut_min_cyc = step_success_tot;
												if (try_mi[j][i] != 0)exp_rec_rate[j][i] = (double)success_mi[j][i] / try_mi[j][i];
												else exp_rec_rate[j][i] = 0;
												success_mi[j][i] = try_mi[j][i] = 0;
												for (k = 0; k < 3; k++)step_try[k] = step_success[k] = 0;
												if (n_mut_here >= step_max)
												{
													if (gen == 0)break;
													else if (stop_li[j][i] == 1)break;
												}
												if (ratio_per_cycle <= ratio_mut_cycle)
												{
													stop_pi[j][i] = 1;
													break;
												}
												fit_mut_prev0 = pop[iter][j][i].fit;
											}
										}
										success_o += success_o_local;
										success_l += success_l_local;
										success_p += success_p_local;
										success_m_tot += success_m_local;
										//if ((gen == 0 && i == 0) || gen > 0)
										{
											fprintf(outlog, "M Pop%d %d %d,%d,%d = %d Try %d M %f H %g I %g F %f", size0[j], i + 1, success_o_local, success_l_local, success_p_local, success_m_local, n_mut_here, pop[iter][j][i].mah, pop[iter][j][i].fpr, pop[iter][j][i].inf, pop[iter][j][i].fit);
											if (gen == 0)fprintf(outlog, "\tOLP %d%d%d RatioP %f", stop_oi[j][i], stop_li[j][i], stop_pi[j][i], ratio[2]);
											fprintf(outlog, "\tGen %d Iter %d\n", gen + 1, iter + 1);
										}
										n_mut_tot += n_mut_here;
										if (success_m_local != 0)success_m++;
									}
								}
								mut_jump[j] = 0;
								for (i = 1; i < MEGE; i++)
								{
									if (pop[iter][j][i].fit >= pop[iter][j][0].fit)mut_jump[j]++;
								}
								qsort((void*)(&pop[iter][j][0]), MEGE, sizeof(pop[iter][j][0]), compare_pop);
							}							
							int trym = 0, sucm = 0;
							for (j = 0; j < CELL; j++)
							{
								for (i = 0; i < MEGE; i++)
								{
									sucm += pop[iter][j][i].sm;
									trym += pop[iter][j][i].tm;
								}
							}
							if (trym != 0)sm0_rate = (double)sucm / trym;
							else sm0_rate = 0;
							if (gen > 0)
							{
								mdo = 0;
								for (j = 0; j < CELL; j++)mdo_cell[j] = 0;
								break;
							}
							else
							{
								int MEGEalf = MEGE / 2;
								double rat0 = 0;
								if (atry[2] > 0)
								{
									rat0 = (double)asuccess[2] / atry[2];
								}
								mdo = 0;
								int m_nostoplo_cell[CELL], m_stopp_cell[CELL];
								for (j = 0; j < CELL; j++)mdo_cell[j] = m_stopp_cell[j] = m_nostoplo_cell[j] = 0;
								for (j = 0; j < CELL; j++)
								{
									if (mut_jump[j] >= MEGEalf)
									{
										mdo = 1;
										mdo_cell[j] = 1;
										fprintf(outlog, "Pop %d Mut jump %d\tSM rate %g\n", size0[j], mut_jump[j], sm0_rate);
									}
									else
									{										
										int m_nostoplo = 0, m_stopp = 0;
										for (i = 0; i < MEGE; i++)
										{
											if ((stop_li[j][i] == 0 || stop_oi[j][i] == 0) && stop_pi[j][i] == 0)
											{
												m_nostoplo++;
											}
											if (stop_pi[j][i] == 1)
											{
												m_stopp++;
											}
										}
										m_nostoplo_cell[j] = m_nostoplo;
										m_stopp_cell[j] = m_stopp;
										if (rat0 == 0)mdo_cell[j] = 0;
										if (m_nostoplo > MEGEalf || rat0 > ratio2_gen0)mdo_cell[j] = 1;
										if (m_stopp >= MEGEalf)mdo_cell[j] = 0;
									}
								}
								fprintf(outlog, "Pop "); for (j = 0; j < CELL; j++)fprintf(outlog, "%d ", size0[j]);
								fprintf(outlog, "\tMut Jump "); for (j = 0; j < CELL; j++)fprintf(outlog, "%d ", mut_jump[j]);
								fprintf(outlog, "\tyLOnP "); for (j = 0; j < CELL; j++)fprintf(outlog, "%d ", m_nostoplo_cell[j]);
								fprintf(outlog, "\tyP "); for (j = 0; j < CELL; j++)fprintf(outlog, " %d ", m_stopp_cell[j]);
								fprintf(outlog, "\tCommon RatioP %f SM rate %f\n", rat0, sm0_rate);
								tnow = time(NULL);
								fprintf(outlog, "Gen %d Iter %d\n", gen + 1, iter + 1);
								fprintf(outlog, "%s", ctime(&tnow));
							}
							mdo = 0;
							for (j = 0; j < CELL; j++)
							{
								if (mdo_cell[j] == 1) { mdo = 1; break; }
							}
						} while (mdo == 1);
						{
							//sm0_rate = 0.001;
						}
					}
					{
						int testar = TestArgv(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], argv[8], argv[9], argv[10], argv[11], argv[12],
							path_fasta, file_for, file_back, reg_max, size_start, size_dif, olen, fp2, octa, path_out, len_peak_max, file_log, outlog);
						if (testar != 0)
						{
							printf("TestArgv After Mut: %s argument %d error\n", argv[testar], testar);
							exit(1);
						}
					}					
										
				//	if (gen == 0)success_m /= (step_max_tot / step_max);
					if (gen <= 1)ratio_thr_r[0] = ratio_thr;		//OL								
					else ratio_thr_r[0] = 1;
					for (j = 0; j < CELL; j++)for (i = 0; i < MEGE; i++)pop[iter][j][i].sr = pop[iter][j][i].tr = 0;					
					for (j = 0; j < CELL; j++)
					{
						int recount = 0;
						double resum = 0;
						for (i = 0; i < MEGE - 1; i++)
						{
							double e1 = Max(ratio_thr, exp_rec_rate[j][i]);
							double e2 = Max(ratio_thr, exp_rec_rate[j][i + 1]);
							resum += success_ri[j][i] * (e1 + e2) / 2;
							recount += success_ri[j][i];
						}
						if (recount != 0)ratio_thr_r[1] = resum / recount;
						else ratio_thr_r[1] = ratio_thr;
					}					
					for (m = 0; m < CELL; m++)fit_after_mut[m] = pop[iter][m][0].fit;
					for (m = 0; m < NREC; m++)success_r1[m] = 0;
					//recombinations																
					fprintf(outlog, "Total mut %d\t", n_mut_tot);
					fprintf(outlog, "SM rate %f", sm0_rate);
					fprintf(outlog, "\t");
					fprintf(outlog, "Mut jump ");
					for (j = 0; j < CELL; j++)fprintf(outlog, " %d", mut_jump[j]);
					fprintf(outlog, "\n");
					fprintf(outlog, "Enter recombinations Gen %d Iteration %d\n", gen+1, iter + 1);
					fprintf(outlog, "\n");
					int loc_rec, loc_rec_tot = 0;					
					int rwei[MEGE];
					double jwei1 = jwei;
					int jmax = MEGE / 2 - 1, jmax1 = jmax - 1;
					rwei[jmax] = 2;
					for (j = jmax1; j >= 0; j--)
					{
						rwei[j] = (int)(rwei[jmax] * jwei1);
						jwei1 *= jwei;
					}
					for (j = jmax; j >= 0; j--)fprintf(outlog, "%d %d\t", j, rwei[j]);
					fprintf(outlog, "\n");
					int pair_all = 0;
					for (j = jmax; j >= 0; j--)
					{
						for (k = j; k < MEGE; k++)
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
					fprintf(outlog, "Total pairs %d Gen %d\tIter %d\n", pair_all, gen + 1,iter+1);					
					{
						int testar = TestArgv(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], argv[8], argv[9], argv[10], argv[11], argv[12],
							path_fasta, file_for, file_back, reg_max, size_start, size_dif, olen, fp2, octa, path_out, len_peak_max, file_log, outlog);
						if (testar != 0)
						{
							printf("TestArgv Before Rec: %s argument %d error\n", argv[testar], testar);
							exit(1);
						}
					}
					for (k = 0; k < pair_all; k++)pair_take[k] = k;
					int n_rec_cycle, sr_step_cycle = sr_step;
					if (gen <= 1) { n_rec_cycle = Min(n_mut_tot / 4, n_rec_cycle_max); }
					else n_rec_cycle = n_rec_cycle_max;
					n_rec_cycle /= pair_all;
					if (n_rec_cycle < 1)n_rec_cycle = 1;
					sr_step_cycle /= pair_all;
					int rec_jump[CELL][MEGE];
					for (m = 0; m < CELL; m++)for (k = 0; k < MEGE; k++)rec_jump[m][k] = 0;					
					for (m = 0; m < CELL; m++)for (k = 0; k < MEGE; k++)fit_pre_rec[m][k] = pop[iter][m][k].fit;
					int sr;
					double ratio_r[2] = { 1,1 };
					int step_rtry[NREC], step_rsuccess[NREC];
					for (k = 0; k < NREC; k++)step_rtry[k] = step_rsuccess[k] = 0;
					int success_r_cycle = 0;
					fprintf(outlog, "Rec %d cycles of %d tries\tTotal %d\tRatioThrOL %.5f RatioThrP %.5f StepR %d\t", n_rec_cycle, pair_all, n_rec_cycle * pair_all, ratio_thr_r[0], ratio_thr_r[1], sr_step_cycle);
					fprintf(outlog, "RE %d", mege_rec);
					fprintf(outlog, "\n");
					
					double fit_rec_prev[CELL], fit_rec_prev0[CELL], mah_rec_prev[CELL], fpr_rec_prev[CELL];
					for (m = 0; m < CELL; m++)fit_rec_prev[m] = pop[iter][m][mege_rec].fit;
					for (m = 0; m < CELL; m++)fit_rec_prev0[m] = pop[iter][m][0].fit;
					for (m = 0; m < CELL; m++)mah_rec_prev[m] = pop[iter][m][0].mah;
					for (m = 0; m < CELL; m++)fpr_rec_prev[m] = pop[iter][m][0].fpr;
					for (m = 0; m < CELL; m++)for (i = 0; i < MEGE; i++)success_ri[m][i] = 0;
					loc_rec = 0;
					for (m = 0; m < CELL; m++)for (i = 0; i < MEGE; i++)pop[iter][m][i].sr = pop[iter][m][i].tr = 0;
					sr = 0;
					int cell1 = CELL - 1;
					for (m = 0; m < pair_all; m++)
					{
						first_cell[m] = second_cell[m] = m % CELL;
					}
					
					//for (sr = 1; sr <= n_rec_cycle; sr++)
					
					do
					{
						sr++;
						BigMixI(pair_take, pair_all);
						BigMixI(second_cell, pair_all);
						for (k = 0; k < pair_all; k++)
						{
							int kk[2];
							for (m = 0; m < 2; m++)kk[m] = pair_d[pair_take[k]][m];
							int pair_size[2];
							{
								pair_size[0] = first_cell[k];
								pair_size[1] = second_cell[k];
							}
							if (pair_size[0] == pair_size[1] && kk[0] == kk[1])continue;
							for (m = 0; m < 2; m++)
							{
								pop[iter][pair_size[m]][kk[m]].get_copy(&det2[m], nseq, reg_max);
							}
							double fit_parent_max;
							fit_parent_max = Max(det2[0].fit, det2[1].fit);
							int nsn[2], nsi[2], num[2];
							int r_cy;
							//if (restart_full == 0)r_cy = rand() % 6;
							//else r_cy = 5;
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
								if (Reco2_Two_dinucleotides_full(&det2[0], &det2[1], n_train[iter], reg_max) == -1)continue;
								else break;
							}
							case(5):
							{
								if (Reco2Peak(&det2[0], &det2[1], n_train[iter], xporti) == -1)continue;
								else break;
							}
							}
							if (r_cy <= 1)
							{
								for (m = 0; m < 2; m++)
								{
									nsn[m] = det2[m].order(nsi[m]);
								}
							}							
							double dd_r[2] = { 0, 0 };
							int gom[2] = { 0, 0 };
							{
								for (m = 0; m < 2; m++)
								{
									for (int t = 0; t < MEGE; t++)
									{
										if (t == pair_d[k][m])continue;
										gom[m] = GomTown(det2[m], pop[iter][pair_size[m]][t], n_train[iter], xporti, check_lpd);
										if (gom[m] == -1)
										{
											break;
										}
									}
									if (gom[m] == -1)break;
								}
							}
							if (gom[0] != 1 && gom[1] != 1)
							{
								for (m = 0; m < 2; m++)
								{
									det2[m].fit = EvalMahFIT(&det2[m], n_train[iter], octa, infc, xporti, seq_real, olen, dav, dcv, octa_prow, pexp, len, outlog);//hoxa_wei, qps, frp, octa_rat,
									if(det2[m].fit>0)pop[iter][pair_size[m]][kk[m]].tr++;
								}
								if (det2[0].fit > 0 || det2[1].fit > 0)step_rtry[r_cy]++;
								{									
									for (m = 0; m < 2; m++)
									{
										if (det2[m].fit > pop[iter][pair_size[m]][kk[m]].fit)
										{
											det2[m].get_copy(&pop[iter][pair_size[m]][kk[m]], nseq, reg_max);
											success_rcell[pair_size[m]]++;
											pop[iter][pair_size[m]][kk[m]].sr++;
											success_r1[r_cy]++;
											step_rsuccess[r_cy]++;
											success_r_cycle++;
											success_ri[pair_size[m]][kk[m]]++;
											loc_rec++;
											if (det2[m].fit > fit_rec_prev[pair_size[m]])
											{
												rec_jump[pair_size[m]][kk[m]]++;
											}
										}
									}									
								}
							}
						}
						int success_r = 0;
						for (j = 0; j < CELL; j++)for (i = 0; i < MEGE; i++)if (success_ri[j][i] > 0)success_r++;
						if (sr % sr_step_cycle == 0)
						{
							int stepr;
							stepr = step_rtry[0] + step_rtry[1] + step_rtry[2] + step_rtry[3] + step_rtry[4];
							if (stepr > 0)ratio_r[0] = (double)(step_rsuccess[0] + step_rsuccess[1] + step_rsuccess[2] + step_rsuccess[3] + step_rsuccess[4]) / stepr;
							else ratio_r[0] = 0;
							stepr = step_rtry[5];
							if (stepr > 0)ratio_r[1] = (double)(step_rsuccess[5]) / stepr;
							else ratio_r[1] = 0;
							for (m = 0; m < CELL; m++)qsort((void*)(&pop[iter][m][0]), MEGE, sizeof(pop[iter][m][0]), compare_pop);
							int tryr = 0, sucr = 0, rec_jump1 = 0;
							for (m = 0; m < CELL; m++)
							{
								for (i = 0; i < MEGE; i++)
								{
									sucr += pop[iter][m][i].sr;
									tryr += pop[iter][m][i].tr;
									if (rec_jump[m][i] > 0)
									{
										rec_jump1++;									
									}
									else break;
								}
							}
							double sr0_rate;
							if (tryr > 0) sr0_rate = (double)sucr / tryr;
							else sr0_rate = 0;
							fprintf(outlog, "Rec %d: %d %d,%d,%d,%d,%d,%d %d OL %f P %f Cell", sr * pair_all, success_r, success_r1[0], success_r1[1], success_r1[2], success_r1[3], success_r1[4], success_r1[5], success_r_cycle, ratio_r[0], ratio_r[1]);
							for (j = 0; j < CELL; j++)fprintf(outlog, " %d", success_rcell[j]);
							fprintf(outlog, "\n");
							for (j = 0; j < CELL; j++)fprintf(outlog, "Pop %d M %f H %g I %f F %f\tRatioSco %f\n", size0[j], pop[iter][j][0].mah, pop[iter][j][0].fpr, pop[iter][j][0].inf, pop[iter][j][0].fit, pop[iter][j][0].fit / fit_rec_prev0[j] - 1);
							for (i = 0; i < CELL; i++)fit_rec_prev[i] = pop[iter][i][mege_rec].fit;
							loc_rec_tot += loc_rec;
							double ratio_per_cycle = 0;
							for (i = 0; i < CELL; i++)
							{
								double rpc = pop[iter][i][0].fit / fit_rec_prev0[i] - 1;
								if (rpc > ratio_per_cycle)ratio_per_cycle = rpc;
							}
							fprintf(outlog, "Gen %d Iteration %d\n", gen+1, iter + 1);
							fprintf(outlog, "L%d SR rate %f Jump %d RatioSco %f\n", loc_rec_tot, sr0_rate, rec_jump1, ratio_per_cycle);
							tnow = time(NULL);
							fprintf(outlog, "%s", ctime(&tnow));
							for (i = 0; i < CELL; i++)for (k = 0; k < MEGE; k++)rec_jump[i][k] = 0;
							if (loc_rec == 0)break;
							if (ratio_r[0] < ratio_thr_r[0] && ratio_r[1] < ratio_thr_r[1])break;
							if (sr0_rate <= sm0_rate)break;
							if (ratio_per_cycle <= ratio_rec_cycle)break;
							if (sr > n_rec_cycle && rec_jump1 <= mege_rec + 1)break;
							for (m = 0; m < NREC; m++)step_rtry[m] = step_rsuccess[m] = 0;
							success_r_cycle = 0;
							loc_rec = 0;
							for (i = 0; i < CELL; i++)fit_rec_prev0[i] = pop[iter][i][0].fit;							
						}
					} 
					while (loc_rec >= 0);										
										
					{
						int testar = TestArgv(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], argv[8], argv[9], argv[10], argv[11], argv[12],
							path_fasta, file_for, file_back, reg_max, size_start, size_dif, olen, fp2, octa, path_out, len_peak_max, file_log, outlog);
						if (testar != 0)
						{
							printf("TestArgv After Rec: %s argument %d error\n", argv[testar], testar);
							exit(1);
						}
					}
					for (i = 0; i < CELL; i++)
					{
						qsort((void*)(&pop[iter][i][0]), MEGE, sizeof(pop[iter][i][0]), compare_pop);
						for (k = 1; k < MEGE; k++)
						{
							if (stop_pi[i][k] == 1 && fit_pre_rec[i][k] > pop[iter][i][k].fit)stop_pi[i][k] = 0;
						}
					}
					double change_level = 0, change_level_rec = 0, change_level_mut = 0;
					for (i = 0; i < CELL; i++)
					{
						double change_level_here = pop[iter][i][0].fit / fit_prev[i] - 1;
						double change_level_rec_here = pop[iter][i][0].fit / fit_after_mut[i] - 1;
						double change_level_mut_here = fit_after_mut[i] / fit_prev[i] - 1;
						if (change_level_here > change_level)change_level = change_level_here;
						if (change_level_rec_here > change_level_rec)change_level_rec = change_level_rec_here;
						if (change_level_mut_here > change_level_mut)change_level_mut = change_level_mut_here;
					}
					double exit_1st = 0.05, exit_2nd = 0.01;
					//double exit_1st = 2, exit_2nd = 0.9;
					if (change_level < exit_2nd && gen > 1)// 3rd iteration at least
					{
						fprintf(outlog, "Go out %d iteration\n", iter + 1);						
						{
							int testar = TestArgv(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], argv[8], argv[9], argv[10], argv[11], argv[12],
								path_fasta, file_for, file_back, reg_max, size_start, size_dif, olen, fp2, octa, path_out, len_peak_max, file_log, outlog);
							if (testar != 0)
							{
								printf("TestArgv Go out %d iteration: %s argument %d error\n", iter+1, argv[testar], testar);
								exit(1);
							}
						}
						big_exit1 = 1;
					}
					else
					{
						fprintf(outlog, "Continue iteration %d\n", iter + 1);
						{
							int testar = TestArgv(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], argv[8], argv[9], argv[10], argv[11], argv[12],
								path_fasta, file_for, file_back, reg_max, size_start, size_dif, olen, fp2, octa, path_out, len_peak_max, file_log, outlog);
							if (testar != 0)
							{
								printf("TestArgv Continue %d iteration: %s argument %d error\n", iter + 1, argv[testar], testar);
								exit(1);
							}
						}
						/*
						if (restart_half > 0 && change_level < exit_1st)
						{
							big_exit1 = 0;
							restart_full++;
							int n_peaks = 20;
							fprintf(outlog, "Restart Full %d\n", restart_full);
							for (i = 0; i < MEGE; i++)for (j = 0; j < CELL; j++)stop_oi[i][j] = stop_li[i][j] = 0;
							int best_alive = MEGE / 3;
							//int best_alive = (int)(MEGE / (3 + 2 * restart_half));//for mege = 30 -> 10,6,4,3,2,2,2
							//if (best_alive < 2)best_alive = 2;
							for (j = 0; j < CELL; j++)
							{
								for (i = best_alive; i < MEGE; i++)
								{
									int best_num = i % best_alive;
									pop[iter][j][best_num].get_copy(&pop[iter][j][i], nseq, reg_max);
									pop[iter][j][i].init_rand_part_hoxa(n_train[iter], xporti, n_peaks, olen, len_octa, len, octa_prowb, octa_prows);
									EvalMahFIT(&pop[iter][j][i], n_train[iter], octa, xporti, seq_real, olen, dav, dcv, octa_prow, len, outlog);//qps, frp, octa_rat,
								}
							}
							for (j = 0; j < CELL; j++)qsort((void*)(&pop[iter][j][0]), MEGE, sizeof(pop[iter][j][0]), compare_pop);
							for (i = 0; i < CELL; i++)for (j = 0; j < MEGE; j++)stop_pi[i][j] = 0;
							for (i = 0; i < CELL; i++)for (j = 0; j < MEGE; j++)stop_oi[i][j] = stop_li[i][j] = 1;
						}
						if (restart_full == 0)
						{
							int best_alive = MEGE / 3;
							//int best_alive = (int)(MEGE / (3 + 2 * restart_half));//for mege = 30 -> 10,6,4,3,2,2,2
							//if (best_alive < 2)best_alive = 2;
							restart_half++;
							fprintf(outlog, "Restart Half %d\n", restart_half);
							int n_peaks = 20;
							for (j = 0; j < CELL; j++)
							{
								for (i = best_alive; i < MEGE; i++)
								{
									int best_num = i % best_alive;
									pop[iter][j][best_num].get_copy(&pop[iter][j][i], nseq, reg_max);
									pop[iter][j][i].init_rand_part_hoxa(n_train[iter], xporti, n_peaks, olen, len_octa, len, octa_prowb, octa_prows);
									EvalMahFIT(&pop[iter][j][i], n_train[iter], octa, xporti, seq_real, olen, dav, dcv, octa_prow, len, outlog);//qps, frp, octa_rat,
									//	stop_pi[j][i] = 0;
								}
							}
							for (j = 0; j < CELL; j++)qsort((void*)(&pop[iter][j][0]), MEGE, sizeof(pop[iter][j][0]), compare_pop);
							for (i = 0; i < CELL; i++)for (j = 0; j < MEGE; j++)stop_pi[i][j] = 0;
							if (restart_half == 1)
							{
								for (i = 0; i < CELL; i++)for (j = 0; j < MEGE; j++)stop_oi[i][j] = stop_li[i][j] = 0;
							}
							else
							{
								for (i = 0; i < CELL; i++)for (j = best_alive; j < MEGE; j++)stop_oi[i][j] = stop_li[i][j] = 0;
							}
						}
						*/
						for (i = 0; i < CELL; i++)for (j = 0; j < MEGE; j++)stop_pi[i][j] = stop_oi[i][j] = stop_li[i][j] = 0;
					}
					gen++;
					for (j = 0; j < CELL; j++)fit_prev[j] = pop[iter][j][0].fit;
					fprintf(outlog, "Gen %d Rat %.5f RatM %.5f RatR %.5f Fit", gen+1, change_level, change_level_mut, change_level_rec);
					for (j = 0; j < CELL; j++)fprintf(outlog, " %.5f", pop[iter][j][0].fit);
					fprintf(outlog, "M %d %d,%d,%d R %d ", success_m, success_o, success_l, success_p, success_r1[4]);
					{
						int sumr = 0;
						for (j = 0; j < CELL; j++)for (i = 0; i < MEGE; i++)sumr += success_ri[j][i];
						if (sumr > 0)
						{
							for (j = 0; j < CELL; j++)
							{
								for (i = 0; i < MEGE; i++)fprintf(outlog, "%.2f ", 100 * (double)success_ri[j][i] / sumr);
								fprintf(outlog, "\n");
							}
						}
					}
					{
						tnow = time(NULL);
						fprintf(outlog, "%s", ctime(&tnow));
					}
					//	if (gen == 1)for (j = 0; j < CELL; j++)for (i = 0; i < MEGE; i++)stop_li[j][i] = stop_oi[j][i] = 1;
					for (j = 0; j < CELL; j++)for (i = 0; i < MEGE; i++)stop_pi[j][i] = 0;
				} while (big_exit1 == 0);
				{
					int testar = TestArgv(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], argv[8], argv[9], argv[10], argv[11], argv[12],
						path_fasta, file_for, file_back, reg_max, size_start, size_dif, olen, fp2, octa, path_out, len_peak_max, file_log, outlog);
					if (testar != 0)
					{
						printf("TestArgv Out iteration %d: %s argument %d error\n", iter + 1, argv[testar], testar);
						exit(1);
					}
				}
				{
					int ic = 0;
					for (k = 0; k < nseq; k++)
					{
						if (xport[k] == 0)xportj[ic++] = k;
					}
				}
				{
					char extmat[6], extmap[6];
					strcpy(extmap,"_cmap");
					strcpy(extmat, "_cmat");
					extmap[5] = '\0', extmat[5] = '\0';
					char file_base[500], name[200];
					memset(name, '\0', sizeof(name));					
					for (i = 0;; i++)
					{
						if (file_for[i] == '.') { name[i] = '\0'; break; }
						if (file_for[i] == '\n') { name[i] = '\0'; break; }
						if (file_for[i] == '\0') { name[i] = '\0'; break; }
						name[i] = file_for[i];
					}
					memset(file_base, '\0', sizeof(file_base));
					strcpy(file_base, path_out);
					strcat(file_base, name);
					strcat(file_base, "_");	
					for (j = 0; j < CELL; j++)
					{			
						for (i = 0; i < MEGE; i++)
						{
							char buf1[5], buf2[5], buf3[5]; 
							memset(buf1, '\0', sizeof(buf1));
							sprintf(buf1, "%d", iter + 1);
							memset(buf2, '\0', sizeof(buf2));
							sprintf(buf2, "%d", size0[j]);
							memset(buf3, '\0', sizeof(buf3));
							sprintf(buf3, "%d", i + 1);							
							/*memset(filemap, '\0', sizeof(filemap));
							strcpy(filemap, file_base);
							strcat(filemap, buf1);
							strcat(filemap, "_");
							strcat(filemap, buf2);
							strcat(filemap, extmap);							
							strcat(filemap, buf3);*/
							memset(filemat, '\0', sizeof(filemat));
							strcpy(filemat, file_base);
							strcat(filemat, buf1);
							strcat(filemat, "_");
							strcat(filemat, buf2);
							strcat(filemat, extmat);
							strcat(filemat, buf3);
							EvalMahControl(&pop[iter][j][i], nseq, nseqb, n_train[iter], n_cntrl[iter], xporti, xportj, seq_real, seq_back, olen, len, lenb, dav, dcv, filemat, outlog, out_roc, out_prc, file_for, fp2,iter,i);//filemap,
						}
					}
				}								
				big_exit1 = 1;
			}	
			{
				int testar = TestArgv(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], argv[8], argv[9], argv[10], argv[11], argv[12],
					path_fasta, file_for, file_back, reg_max, size_start, size_dif, olen, fp2, octa, path_out, len_peak_max, file_log, outlog);
				if (testar != 0)
				{
					printf("TestArgv Out EvalMahControl %d: %s argument %d error\n", iter + 1, argv[testar], testar);
					exit(1);
				}
			}
			for (i = 0; i < MEGE; i++)
			{
				for (j = 0; j < CELL; j++)
				{
					fprintf(out_auc, "%s\t%d\t%d\t%d", file_for, olen, pop[0][j][i].size, i + 1);
					for (k = 0; k < 2; k++)fprintf(out_auc, "\t%f", pop[k][j][i].aucroc);
					for (k = 0; k < 2; k++)fprintf(out_auc, "\t%f", pop[k][j][i].auprc);
					for (k = 0; k < 2; k++)fprintf(out_auc, "\t%f", pop[k][j][i].mah);
					for (k = 0; k < 2; k++)fprintf(out_auc, "\t%f", pop[k][j][i].fpr);
					for (k = 0; k < 2; k++)fprintf(out_auc, "\t%f", pop[k][j][i].inf);
					fprintf(out_auc, "\n");
				}
			}
			fclose(out_prc);
			fclose(out_roc);
			fclose(out_auc);
		}
	}
	fclose(outlog);
	for (iter = 0; iter < iteration; iter++)
	{
		for (j = 0; j < CELL; j++)for (i = 0; i < MEGE; i++)pop[iter][j][i].mem_out();
		delete[]pop[iter];
	}
	delete[] pop;
	for (i = 0; i < 2; i++)det2[i].mem_out();
	det1.mem_out();
	delete[] pair_take;
	delete[] first_cell;
	delete[] second_cell;
	for (i = 0; i < pair_all_max; i++)
	{
		delete[] pair_d[i];
	}
	delete[] pair_d;
	for (k = 0; k < 2; k++)
	{
		for (i = 0; i < nseq; i++)
		{
			delete[] seq_back[k][i];
		}
		delete[] seq_back[k];
	}
	delete[] seq_back;
	for (k = 0; k < 2; k++)
	{
		for (i = 0; i < nseq; i++)
		{
			delete[] seq_real[k][i];
		}
		delete[] seq_real[k];
	}
	delete[] seq_real;
	for (i = 0; i < reg_max; i++)
	{
		delete[] dcv[i];
	}
	delete[] dcv;
	for (i = 0; i < reg_max; i++)
	{
		delete[] dav[i];
	}
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
	for (i = 0; i < nseq; i++)
	{
		delete[] octa_prows[i];
	}
	delete[] octa_prows;
	delete[] octa_pro1p;
	delete[] len;
	delete[] lenb;
	delete[] len_octa;
	delete[] thr_octa;
	delete[] xport;
	delete[] xporti;
	delete[] xportj;
	delete[] n_train;
	delete[] n_cntrl;
	delete[] octa_rat;
	return 0;
}