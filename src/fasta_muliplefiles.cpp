#define _CRT_SECURE_NO_WARNINGS

#include  <stdio.h>
#include  <stdlib.h>
#include  <string.h>
#include  <math.h>
#include  <time.h>

#define SEQLEN 5000
//#define DELTA 5000 // perekryvanie mejdu posl-tyami

void DelHole(char *str)
{
   char *hole;

   hole=strstr(str,"\n");
   if(hole!=NULL) *hole=0;
}

int main(int argc, char *argv[])
{
int nseq=0;
FILE *in, *out;
char head[SEQLEN], basename[80], file[80], str[50];// , d[DELTA];
int string_len=50000;
int shift=200;
char ext[]=".fa";

if(argc!=4)
{
	printf("Command line error: %s 1 file_input_fasta 2file_out_name 3int output filenames: 1 = 1,2,3, etc; 0 = from headers of sequences", argv[0]);
    exit(1);
}
if((in=fopen(argv[1],"rt"))==NULL)
{
	printf("Input file %s can't be opened", argv[1]);
	exit(1);
}
int length=0;
strcpy(basename,argv[2]);
int mode = atoi(argv[3]);
char c, symbol;
symbol=getc(in);
rewind(in);
int part;
out = NULL;
char test1[10], test2[10], sep1 = ' ', sep2 = '\t';
strcpy(test1, "chr");
strcpy(test2, "Chr");
  while((c=getc(in))!=EOF)
  {
	if(c==symbol)
	{
		nseq++;
		length=0;
		part=1;
		fgets(head,sizeof(head),in);
		DelHole(head);
		if(nseq>1)fclose(out);
		memset(file,0,sizeof(file));
		strcpy(file,basename);				
		memset(str, 0, sizeof(str));
		if (mode == 1)
		{
			sprintf(str, "%d", nseq);
		}
		else
		{			
			int slen = 0;
			int hlen = strlen(head);
			int j;			
			for (j = 0; j < hlen; j++)
			{
				if (head[j] == sep1 || head[j] == sep2)
				{					
					continue;
				}
				break;
			}	
			int j0 = j;
			//if (hlen - j0 >= 3)
			{
				if (strncmp(&head[j0], test1, 3) == 0 || strncmp(&head[j0], test2, 3) == 0)j0 += 3;
			}
			for (j = j0;j<hlen; j++)
			{
				if (head[j] == sep1 || head[j] == sep2)
				{
					if (slen == 0)
					{
						printf("Error header recognition file %s header %s\n", argv[1], head);
						exit(1);
					}
					break;
				}
				else
				{
					str[slen++] = head[j];
				}
			}
			if (slen == 0)
			{
				printf("Error header recognition file %s header %s\n", argv[1], head);
				exit(1);
			}
			str[slen] = '\0';
		}
		strcat(file, str);
		strcat(file,ext);
		if((out=fopen(file,"w+t"))==NULL)
		{
			puts("Output file can't be opened");
			exit(1);
		}
		//if(nseq!=1)fprintf(out,"\n");	
		fprintf(out,">%s\n",head);	
		continue;
	}
//	if(c=='\n')continue;	
	length++;	
	fprintf(out,"%c",c);
	/*if(length==string_len)
	{
		part++;
		memset(d,0,sizeof(d));
		long pp=ftell(out);
		fseek(out,pp-(long)shift,SEEK_SET);
		for(int i=0;i<shift;i++)
		{
			c2=getc(out);
			d[i]=c2;
		}		
		fseek(out,pp,SEEK_SET);
		d[shift]='\0';
		fprintf(out,"\n>%s\tpart%d\n%s",head,part,d);
		length=0;
	}*/
  }
fclose(in);
//fprintf(out,"\n");
fclose(out);
return 1;
}