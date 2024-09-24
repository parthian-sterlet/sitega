#define _CRT_SECURE_NO_WARNINGS

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
//#include  <conio.h>
#include "chromosome_name.h"
#define SEQLEN 20001
#define NCHR 100

void DelHole(char* str)
{
	char* hole;

	hole = strstr(str, "\n");
	if (hole != NULL) *hole = 0;
}
char TransSym(char c)
{
	int cc = int(c);
	if (cc >= 97) c = char(cc - 32);
	return(c);
}

int main(int argc, char* argv[])
{
	FILE* in, * out;
	char path_fasta[500], genome[10], d[SEQLEN];
	char filesta[500], fileoend[10], fileiend[10], filei[500], fileo[500];

	if (argc != 3)
	{
		printf("Command line error: %s 1path_genome  2char genome (hg38 mm10 rn6 zf11 dm6 ce235 sc64 sch294 at10 gm21 zm73 mp61)", argv[0]);
		exit(1);
	}
	strcpy(path_fasta, argv[1]);
	strcpy(genome, argv[2]);
	
	int nc, i, n_chr;
	char name_chr[NCHR][10];	
	int genome_rec = 0;
	if (strcmp(genome, "at10") == 0)
	{
		genome_rec = 1;
		n_chr = n_chr_at;
		for (i = 0; i < n_chr; i++)
		{
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
	strcpy(filesta, path_fasta);
	strcat(filesta, "chr");
	strcpy(fileiend, ".fa");
	strcpy(fileoend, ".plain");

	int len, fl, n;
	for (nc = 0; nc < n_chr; nc++)
	{
		strcpy(filei, filesta);		
		strcat(filei, name_chr[nc]);
		strcpy(fileo, filei);
		strcat(filei, fileiend);
		strcat(fileo, fileoend);
		if ((in = fopen(filei, "rt")) == NULL)
		{
			printf("Input file %s can't be opened!",filei);
			exit(1);
		}
		if ((out = fopen(fileo, "wt")) == NULL)
		{
			printf("Output file %s can't be opened!", fileo);
			exit(1);
		}
		char c, symbol = '>';
		char alfavit4[] = "atgcATGC";
		char alfavit11[] = "wrmkysbvhdnWRMKYSBVHDN";
		char alfavit15maska[] = "atgcwrmkysbvhdnxATGCWRMKYSBVHDNX";
		char alfavit15maska_low[] = "atgcwrmkysbvhdnx";
		char nka = 'N';
		n = 0;
		fl = 1;
		len = 0;
		n = 0;
		do
		{
			c = getc(in);
			if (c == EOF)fl = -1;
			if (c == symbol || fl == -1)
			{
				if (n > 0)
				{
					fprintf(out, "\n");
				}
				if (c == symbol)n++;
				if (fl == 1)
				{
					fgets(d, sizeof(d), in);
					continue;
				}
			}
			if (c == '\n')
			{
				continue;
			}
			if (c == '\t')continue;
			if (c == ' ')continue;
			if (c == '\r')continue;
			if (fl == 1)
			{
				{
					if (strchr(alfavit4, c) != NULL)
					{
						c = TransSym(c);//segal
						fprintf(out, "%c", c);
						len++;
					}
					else
					{
						if (strchr(alfavit11, c) != NULL)
						{
							//	fprintf(out,"%c",c);					
							fprintf(out, "%c", nka);//segal
							//	printf("%c",c);
							len++;
						}
						else
						{
							{
								printf("Wrong %c symbol found!", c);
								exit(1);
							}
						}
					}
					if (len % 10000000 == 0)printf("\b\b\b\b\b\b\b\b\b\b%d", len);
				}
			}
		} while (fl == 1);
		printf("\b\b\b\b\b\b\b\b\b\b");
		printf("%s %d bp completed\n", name_chr[nc], len);
		fclose(in);
		fclose(out);
	}
	return 1;
}