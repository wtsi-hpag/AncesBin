/****************************************************************************
 ****************************************************************************
 *                                                                          *
 *  Copyright (C) 2017  Genome Research Ltd.                                *
 *                                                                          *
 *  Author: Zemin Ning (zn1@sanger.ac.uk)                                   *
 *                                                                          *
 *  This file is part of bin10x pipeline.                                 *
 *                                                                          *
 *  Scaff10x is a free software: you can redistribute it and/or modify it   *
 *  under the terms of the GNU General Public License as published by the   *
 *  Free Software Foundation, either version 3 of the License, or (at your  *
 *  option) any later version.                                              *
 *                                                                          *
 *  This program is distributed in the hope that it will be useful, but     *
 *  WITHOUT ANY WARRANTY; without even the implied warranty of              *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU        *
 *  General Public License for more details.                                *
 *                                                                          *
 *  You should have received a copy of the GNU General Public License along *
 *  with this program.  If not, see <http://www.gnu.org/licenses/>.         *
 *                                                                          *
 ****************************************************************************
 ****************************************************************************/
/****************************************************************************/

 
#include <math.h>
#include <values.h>
#include <stdio.h>
#include <netinet/in.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include "fasta.h"

#define GT '>'
#define GT4 (((((GT<<8)+GT)<<8)+GT)<<8)+GT

#define ENDS_EXTRA 0
#define PADCHAR '-'
#define MAX_N_BRG 50000 
#define MAX_N_ROW 40000 
#define Max_N_NameBase 400 
#define Max_N_NameBase2 400 
#define Max_N_Pair 100
static char bindir[2000];
static char tmpdir[2000];
static char **S_Name,**R_Name,**R1_Name,**R2_Name,**R1_Name1,**R2_Name1,**R1_Name2,**R2_Name2;
static int *reads_tag1,*reads_tag2,*insert_dev,*core_list,*ctg_list,*ctg_head,*read2contig;
static int *readIndex;

/* SSAS default parameters   */
static int n_group=0;
static int num_reads=10;
static int num1_reads=8;
static int num2_reads=10;
static int num_links=8;
static int num1_links=8;
static int num2_links=10;
static int len_block=50000;
static int len_edges=50000;
static int len_matrx=2000;
static int min_edge = 5;
static int file_tag = 0;
static int data_dir = 0;
static int gap_len = 100;
static int uplinks = 50;
static int mscore = 20;
static int n_longread = 1;
static int run_align = 1;
static int plot_10x = 0;
static int min_len = 3000;
static int sam_flag = 0;

typedef struct
{
       int foffset;
       int fsindex;
} SIO;

void
RunSystemCommand(char *cmd)
{
    int ret;
    if ((ret = system(cmd)) != 0)
    {
        fprintf(stderr, "Error running command: %s\n", cmd);
        exit(EXIT_FAILURE);
    }
}


int main(int argc, char **argv)
{
    int i,j,nSeq,args,nseq;
    int n_files,n_r1,n_r2;
    char *st,*ed;
    char **cmatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch);
    void ArraySort_String2(int n,char **Pair_Name,int *brr);
    fasta *seq;
    FILE *fp,*namef,*namef2;
    int size_file;
    int m_score,n_nodes,n_debug,num_sigma;
    void decodeReadpair(int nSeq);
    void HashFasta_Head(int i, int nSeq);
    void HashFasta_Table(int i, int nSeq);
    void Search_SM(fasta *seq,int nSeq);
    void Assemble_SM(int arr,int brr);
    void Read_Index(int nSeq,char *namefile);
    void Read_Group(fasta *seq,int nSeq,int nRead,int cindex);
    void File_Output(int aaa);
    void Memory_Allocate(int arr);
    char tempa[2000],tempc[2000],syscmd[2000],workdir[2000],line[2000],readname[2000],tagname[10];
    char file_bamseq[2000],file_matR1[2000],file_matR2[2000],file_patR1[2000],file_patR2[2000];
    char file_reads[2000],file_seqdir[2000],samname[500],bamname[500],toolname[500],datname[500];
    int systemRet = system (syscmd);
    int systemChd = chdir(tmpdir);
    pid_t pid;

    seq=NULL;
    
    if(argc < 2)
    {
         printf("Program: bin10X - bin 10X reads based on two haplotype assemblies\n");
         printf("Version: 1.0\n");
         printf("\n");
         
         printf("Usage: %s -nodes 30 <input_sorted_bam> <Input_data_file> <Output_director_binnedreads>\n",argv[0]);
         printf("       nodes    (30)   - number of CPUs requested\n");
         printf("       Input_sorted_bam: read name sorted longranger BAM file\n");
         printf("       Input_data_file:  a text file to point the locations of the reads in paired files\n");
         printf("       Output_director:  a director contained all the binned reads\n");
         exit(1);
    }

    m_score = 50;
    n_nodes = 20;
    n_debug = 1;

    strcpy(toolname,"bwa");
    nSeq=0;
    args=1;
    for(i=1;i<argc;i++)
    {
       if(!strcmp(argv[i],"-gap"))
       {
         sscanf(argv[++i],"%d",&gap_len); 
         args=args+2;
       }
       else if(!strcmp(argv[i],"-nodes"))
       {
         sscanf(argv[++i],"%d",&n_nodes);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-data"))
       {
         data_dir = 1;
//         run_align = 1;
         sscanf(argv[++i],"%s",datname);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-align"))
       {
         memset(toolname,'\0',500);
         sscanf(argv[++i],"%s",toolname);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-score"))
       {
         sscanf(argv[++i],"%d",&mscore);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-help"))
       {
         printf("Usage: %s -nodes 30 <input_sorted_bam> <Input_data_file> <Output_director_binnedreads>\n",argv[0]);
         printf("       nodes    (30)   - number of CPUs requested\n");
         printf("       Input_sorted_bam: read name sorted longranger BAM file\n");
         printf("       Input_data_file:  a text file to point the locations of the reads in paired files\n");
         printf("       Output_director:  a director contained all the binned reads\n");
         exit(1);
       }
       else if(!strcmp(argv[i],"-debug"))
       {
         sscanf(argv[++i],"%d",&n_debug);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-file"))
       {
         sscanf(argv[++i],"%d",&file_tag); 
         args=args+2;
       }
    }

    pid = getpid();
    memset(tempa,'\0',2000);
    if (!getcwd(tempa, sizeof(tempa)))
    {
      exit(1);
    } 
    memset(tmpdir,'\0',2000);
    memset(workdir,'\0',2000);
    sprintf(tmpdir,"%s/",tempa);
    sprintf(workdir,"%s/tmp_rununik_%d/",tempa,pid);

    memset(syscmd,'\0',2000);
    sprintf(syscmd,"mkdir %s",workdir);
    RunSystemCommand(syscmd);

    if(chdir(workdir) == -1)
    {
      printf("System command error: chdir\n");
      exit(EXIT_FAILURE);
    }
     
    st = argv[0];
    ed = strrchr(argv[0],'/');
    memset(tempc,'\0',2000);
    strncpy(tempc,argv[0],ed-st);
    memset(bindir,'\0',2000);
    sprintf(bindir,"%s/seqbin-bin",tempc);

    memset(file_bamseq,'\0',2000);
    memset(file_reads,'\0',2000);
    memset(file_seqdir,'\0',2000);
    memset(file_matR1,'\0',2000);
    memset(file_matR2,'\0',2000);
    memset(file_patR1,'\0',2000);
    memset(file_patR2,'\0',2000);

    sprintf(file_bamseq,"%s/%s",tempa,argv[args]);
    sprintf(file_reads,"%s/%s",tempa,argv[args+1]);
    sprintf(file_seqdir,"%s/%s",tempa,argv[args+2]);
    sprintf(file_matR1,"%s/%s_1.MAT.gz",tempa,argv[args+2]);
    sprintf(file_matR2,"%s/%s_2.MAT.gz",tempa,argv[args+2]);
    sprintf(file_patR1,"%s/%s_1.PAT.gz",tempa,argv[args+2]);
    sprintf(file_patR2,"%s/%s_2.PAT.gz",tempa,argv[args+2]);

    memset(syscmd,'\0',2000);
    sprintf(syscmd,"mkdir %s",file_seqdir);
    RunSystemCommand(syscmd);

    if((namef = fopen(file_bamseq,"r")) == NULL)
    {
      printf("File not in the working directory!\n");
      if((namef = fopen(argv[args],"r")) == NULL)
      {
        printf("File %s not found and please copy it to your working directory!\n",argv[args]);
        exit(1);
      }
      else
      {
        memset(file_bamseq,'\0',2000);
        strcpy(file_bamseq,argv[args]);
        printf("Input sorted bam file1: %s\n",file_bamseq);
      }
    }
    else
    {
      printf("Input sorted bam file2: %s\n",file_bamseq);
    } 

    memset(syscmd,'\0',2000);
//    sprintf(syscmd,"%s/samtools view %s | awk '%s' > align.dat",bindir,file_bamseq,"($2<100){print $1,$2,$3,$5}",bindir);
    sprintf(syscmd,"%s/samtools view %s | awk '%s' | %s/read-binning.pl > align.dat",bindir,file_bamseq,"($2<200){print $1,$2,$3,$5}",bindir);
    RunSystemCommand(syscmd);
    
    memset(syscmd,'\0',2000);
    sprintf(syscmd,"egrep MAT align.dat | awk '{print $1}' > name-MAT.dat");
    RunSystemCommand(syscmd);
    
    memset(syscmd,'\0',2000);
    sprintf(syscmd,"egrep PAT align.dat | awk '{print $1}' > name-PAT.dat");
    RunSystemCommand(syscmd);

    if((namef = fopen(file_reads,"r")) == NULL)
    {
      printf("File not in the working directory!\n");
      if((namef = fopen(argv[args+1],"r")) == NULL)
      {
        printf("File %s not found and please copy it to your working directory!\n",argv[args+1]);
        exit(1);
      }
      else
      {
        memset(file_reads,'\0',2000);
        strcpy(file_reads,argv[args+1]);
        printf("Input reads file: %s\n",file_reads);
      }
    }
    else
    {
      printf("Input reads file: %s\n",file_reads);
    }

    n_files = 0;
    while(!feof(namef))
    {
      fgets(line,2000,namef);
      if(feof(namef)) break;
      n_files++;
    }
    fclose(namef);

    if((reads_tag1 = (int *)calloc(n_files,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - reads_tag\n");
      exit(1);
    }
    if((reads_tag2 = (int *)calloc(n_files,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - reads_tag\n");
      exit(1);
    }

    R1_Name=cmatrix(0,n_files+10,0,2000);
    R2_Name=cmatrix(0,n_files+10,0,2000);
    R1_Name1=cmatrix(0,n_files+10,0,2000);
    R2_Name1=cmatrix(0,n_files+10,0,2000);
    R1_Name2=cmatrix(0,n_files+10,0,2000);
    R2_Name2=cmatrix(0,n_files+10,0,2000);
    S_Name=cmatrix(0,1000,0,Max_N_NameBase);
    R_Name=cmatrix(0,1000,0,Max_N_NameBase);

    if((namef = fopen(file_reads,"r")) == NULL)
    {
      printf("ERROR main:: reads group file \n");
      exit(1);
    }
    i = 0;

    n_r1 = 0;
    n_r2 = 0;
    while(fscanf(namef,"%s",readname)!=EOF)
    {
        int clen = strlen(readname);
        st = strchr(readname,'=');
        ed = strrchr(readname,'/');
        if((strncmp(readname,"q1=",3))==0)
        {
          strcpy(R1_Name[n_r1],st+1);
          strcpy(R1_Name1[n_r1],ed+1);
          ed = strrchr(readname,'.');
          memset(tagname,'\0',10);
          if(ed != NULL)
          {
            strcpy(tagname,ed);
            if((strncmp(tagname,".gz",3))==0)
              reads_tag1[n_r1] = 1;
            else
              reads_tag1[n_r1] = 2;
          }
          else
          {
              reads_tag1[n_r1] = 2;
          }
      printf("file: %d %d %s %s\n",clen,n_r1,R1_Name1[n_r1],readname);
          n_r1++;
        }
        else if((strncmp(readname,"q2=",3))==0) 
        {
          strcpy(R2_Name[n_r2],st+1);
          strcpy(R2_Name1[n_r2],ed+1);
          ed = strrchr(readname,'.');
          memset(tagname,'\0',10);
          if(ed != NULL)
          {
            strcpy(tagname,ed);
            if((strncmp(tagname,".gz",3))==0)
              reads_tag2[n_r2] = 1;
            else
              reads_tag2[n_r2] = 2;
          }
          else
          {
              reads_tag2[n_r2] = 2;
          }
//      printf("file: %d %s \n",n_r2,R2_Name[n_r2]);
          n_r2++;
        }
        i++;
    }
    fclose(namef);
    if(n_r1 != n_r2)
    {
      printf("Number of reads files: %d NOT equal to read2 files: %d\n",n_r1,n_r2);
      memset(syscmd,'\0',2000);
      sprintf(syscmd,"rm -rf %s > /dev/null",workdir);
      RunSystemCommand(syscmd);
      exit(1);
    }
   
    for(i=0;i<n_r1;i++)
    {

       memset(syscmd,'\0',2000);
       sprintf(syscmd,"cp %s %s",R1_Name[i],workdir);
       RunSystemCommand(syscmd);
       if(reads_tag1[i]==1)
       {
         int clen = strlen(R1_Name1[i]);
         memset(syscmd,'\0',2000);
         sprintf(syscmd,"%s/pigz -dc %s > reads_10X_1.fq ",bindir,R1_Name1[i]);
         RunSystemCommand(syscmd);
         strncpy(R1_Name2[i],R1_Name1[i],clen-3);
       }
       else
         strcpy(R1_Name2[i],R1_Name1[i]);

       memset(syscmd,'\0',2000);
       sprintf(syscmd,"%s/seqbin_split -split 50000000 split-reads_1 reads_10X_1.fq > try.out",bindir);
       RunSystemCommand(syscmd);

       memset(syscmd,'\0',2000);
       sprintf(syscmd,"ls -lrt split-reads_1_* | awk '{print $9}' > name1.dat");
       RunSystemCommand(syscmd);

       nseq=0;
       if((namef = fopen("name1.dat","r")) == NULL)
       {
         printf("ERROR main:: args \n");
         exit(1);
       }
       while(!feof(namef))
       {
         fgets(line,2000,namef);
         if(feof(namef)) break;
         nseq++;
       }
       fclose(namef);    

       if((namef = fopen("name1.dat","r")) == NULL)
       {
         printf("ERROR main:: reads group file \n");
         exit(1);
       }

       for(j=0;j<nseq;j++)
          memset(S_Name[j],'\0',Max_N_NameBase);
       
/*  read the alignment files         */
       j=0;
       while(fscanf(namef,"%s",S_Name[j])!=EOF)
       {
          j++;
       }
       fclose(namef);
     
       for(j=0;j<nseq;j++)
       {
          memset(syscmd,'\0',2000);
          sprintf(syscmd,"%s/seqbin_reads name-MAT.dat %s %s.MAT > try.out",bindir,S_Name[j],S_Name[j]);
          RunSystemCommand(syscmd);

          memset(syscmd,'\0',2000);
          sprintf(syscmd,"%s/seqbin_reads name-PAT.dat %s %s.PAT > try.out",bindir,S_Name[j],S_Name[j]);
          RunSystemCommand(syscmd);
       }
       
       memset(syscmd,'\0',2000);
       sprintf(syscmd,"cat split-reads_1_0*.MAT > MAT-%s",R1_Name2[i]);
       RunSystemCommand(syscmd);
 
       memset(syscmd,'\0',2000);
       sprintf(syscmd,"%s/seqbin_bcodes MAT-%s MAT-%s.I1 > try.out",bindir,R1_Name2[i],R1_Name2[i]);
       RunSystemCommand(syscmd);

       memset(syscmd,'\0',2000);
       sprintf(syscmd,"%s/pigz -p %d MAT-%s > try.out",bindir,n_nodes,R1_Name2[i]);
       RunSystemCommand(syscmd);

       memset(syscmd,'\0',2000);
       sprintf(syscmd,"%s/pigz -p %d MAT-%s.I1 > try.out",bindir,n_nodes,R1_Name2[i]);
       RunSystemCommand(syscmd);

       memset(syscmd,'\0',2000);
       sprintf(syscmd,"cat split-reads_1_0*.PAT > PAT-%s",R1_Name2[i]);
       RunSystemCommand(syscmd);
 
       memset(syscmd,'\0',2000);
       sprintf(syscmd,"%s/seqbin_bcodes PAT-%s PAT-%s.I1 > try.out",bindir,R1_Name2[i],R1_Name2[i]);
       RunSystemCommand(syscmd);

       memset(syscmd,'\0',2000);
       sprintf(syscmd,"%s/pigz -p %d PAT-%s > try.out",bindir,n_nodes,R1_Name2[i]);
       RunSystemCommand(syscmd);

       memset(syscmd,'\0',2000);
       sprintf(syscmd,"%s/pigz -p %d PAT-%s.I1 > try.out",bindir,n_nodes,R1_Name2[i]);
       RunSystemCommand(syscmd);

       memset(syscmd,'\0',2000);
       sprintf(syscmd,"mv MAT-%s.gz %s ",R1_Name2[i],file_seqdir);
       RunSystemCommand(syscmd);
 
       memset(syscmd,'\0',2000);
       sprintf(syscmd,"mv MAT-%s.I1.gz %s ",R1_Name2[i],file_seqdir);
       RunSystemCommand(syscmd);
 
       memset(syscmd,'\0',2000);
       sprintf(syscmd,"mv PAT-%s.gz %s ",R1_Name2[i],file_seqdir);
       RunSystemCommand(syscmd);
 
       memset(syscmd,'\0',2000);
       sprintf(syscmd,"mv PAT-%s.I1.gz %s ",R1_Name2[i],file_seqdir);
       RunSystemCommand(syscmd);
 
       memset(syscmd,'\0',2000);
       sprintf(syscmd,"rm -rf reads_10X_1.fq split-reads_1_0* %s > try.out",R1_Name1[i]);
       RunSystemCommand(syscmd);
    }
 
    for(i=0;i<n_r2;i++)
    {

       memset(syscmd,'\0',2000);
       sprintf(syscmd,"cp %s %s",R2_Name[i],workdir);
       RunSystemCommand(syscmd);
       if(reads_tag2[i]==1)
       {
         int clen = strlen(R2_Name1[i]);
         memset(syscmd,'\0',2000);
         sprintf(syscmd,"%s/pigz -dc %s > reads_10X_2.fq ",bindir,R2_Name1[i]);
         RunSystemCommand(syscmd);
         strncpy(R2_Name2[i],R2_Name1[i],clen-3);
       }
       else
         strcpy(R2_Name2[i],R2_Name1[i]);

       memset(syscmd,'\0',2000);
       sprintf(syscmd,"%s/seqbin_split -split 50000000 split-reads_2 reads_10X_2.fq > try.out",bindir);
       RunSystemCommand(syscmd);

       memset(syscmd,'\0',2000);
       sprintf(syscmd,"ls -lrt split-reads_2_* | awk '{print $9}' > name2.dat");
       RunSystemCommand(syscmd);

       nseq=0;
       if((namef = fopen("name2.dat","r")) == NULL)
       {
         printf("ERROR main:: args \n");
         exit(1);
       }
       while(!feof(namef))
       {
         fgets(line,2000,namef);
         if(feof(namef)) break;
         nseq++;
       }
       fclose(namef);    
       R_Name=cmatrix(0,nseq,0,Max_N_NameBase);

       if((namef = fopen("name2.dat","r")) == NULL)
       {
         printf("ERROR main:: reads group file \n");
         exit(1);
       }

       for(j=0;j<nseq;j++)
          memset(R_Name[j],'\0',Max_N_NameBase);

/*  read the alignment files         */
       j=0;
       while(fscanf(namef,"%s",R_Name[j])!=EOF)
       {
          j++;
       }
       fclose(namef);
     
       for(j=0;j<nseq;j++)
       {
          memset(syscmd,'\0',2000);
          sprintf(syscmd,"%s/seqbin_reads name-MAT.dat %s %s.MAT > try.out",bindir,R_Name[j],R_Name[j]);
          RunSystemCommand(syscmd);

          memset(syscmd,'\0',2000);
          sprintf(syscmd,"%s/seqbin_reads name-PAT.dat %s %s.PAT > try.out",bindir,R_Name[j],R_Name[j]);
          RunSystemCommand(syscmd);
       }
       
       memset(syscmd,'\0',2000);
       sprintf(syscmd,"cat split-reads_2_0*.MAT > MAT-%s",R2_Name2[i]);
       RunSystemCommand(syscmd);
 
       memset(syscmd,'\0',2000);
       sprintf(syscmd,"%s/pigz -p %d MAT-%s > try.out",bindir,n_nodes,R2_Name2[i]);
       RunSystemCommand(syscmd);

       memset(syscmd,'\0',2000);
       sprintf(syscmd,"cat split-reads_2_0*.PAT > PAT-%s",R2_Name2[i]);
       RunSystemCommand(syscmd);
 
       memset(syscmd,'\0',2000);
       sprintf(syscmd,"%s/pigz -p %d PAT-%s > try.out",bindir,n_nodes,R2_Name2[i]);
       RunSystemCommand(syscmd);

       memset(syscmd,'\0',2000);
       sprintf(syscmd,"mv MAT-%s.gz %s ",R2_Name2[i],file_seqdir);
       RunSystemCommand(syscmd);
 
       memset(syscmd,'\0',2000);
       sprintf(syscmd,"mv PAT-%s.gz %s ",R2_Name2[i],file_seqdir);
       RunSystemCommand(syscmd);
 
       memset(syscmd,'\0',2000);
       sprintf(syscmd,"rm -rf reads_10X_2.fq split-reads_2_0* %s > try.out",R2_Name1[i]);
       RunSystemCommand(syscmd);
    }

    if(n_debug == 0)
    {
      memset(syscmd,'\0',2000);
      sprintf(syscmd,"rm -rf * > /dev/null");
      RunSystemCommand(syscmd);
      chdir(tmpdir);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"rm -rf %s > /dev/null",workdir);
      RunSystemCommand(syscmd);
    }
    return EXIT_SUCCESS;

}
/* end of the main */



#define SWAP(a,b) temp=(a);(a)=b;(b)=temp;

/*   Subroutine to sort an array arr[0,...,n-1] into ascending order while
     making the corresponding reaarangement of the array brr[0,...,n-1]
     by the use of Quicksort (Sedgwick, R. 1978, Communications o fthe ACM,
     vol. 21, pp. 847-857) also see Numerical Recipes in C                  */  

/* =============================== */
void ArraySort_Long(int n, B64_long *arr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,NSTACK=50,istack[NSTACK];
     B64_long a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]<=a) break;
                arr[i+1]=arr[i];
             }
             arr[i+1]=a;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);

          if(arr[m]>arr[ir])
          {
            SWAP(arr[m],arr[ir]);
          }

          if(arr[m+1]>arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
          }

          if(arr[m]>arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          for(;;)
          {
             do i++; while (arr[i]<a);
             do j--; while (arr[j]>a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}


/* =============================== */
void ArraySort_Int(int n, int *arr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,NSTACK=50,istack[NSTACK];
     int a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]<=a) break;
                arr[i+1]=arr[i];
             }
             arr[i+1]=a;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);

          if(arr[m]>arr[ir])
          {
            SWAP(arr[m],arr[ir]);
          }

          if(arr[m+1]>arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
          }

          if(arr[m]>arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          for(;;)
          {
             do i++; while (arr[i]<a);
             do j--; while (arr[j]>a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}


/* =============================== */
void ArraySort_Mix(int n, B64_long *arr, int *brr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,NSTACK=50,istack[NSTACK];
     B64_long a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             b=brr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]<=a) break;
                arr[i+1]=arr[i];
                brr[i+1]=brr[i];
             }
             arr[i+1]=a;
             brr[i+1]=b;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);
          SWAP(brr[k],brr[m+1]);

          if(arr[m]>arr[ir])
          {
            SWAP(arr[m],arr[ir]);
            SWAP(brr[m],brr[ir]);
          }

          if(arr[m+1]>arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
            SWAP(brr[m+1],brr[ir]);
          }

          if(arr[m]>arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
            SWAP(brr[m],brr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          b=brr[m+1];
          for(;;)
          {
             do i++; while (arr[i]<a);
             do j--; while (arr[j]>a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
             SWAP(brr[i],brr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
          brr[m+1]=brr[j];
          brr[j]=b;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}

/* =============================== */
void ArraySort_Int2(int n, int *arr, int *brr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,NSTACK=50,istack[NSTACK];
     int a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             b=brr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]<=a) break;
                arr[i+1]=arr[i];
                brr[i+1]=brr[i];
             }
             arr[i+1]=a;
             brr[i+1]=b;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);
          SWAP(brr[k],brr[m+1]);

          if(arr[m]>arr[ir])
          {
            SWAP(arr[m],arr[ir]);
            SWAP(brr[m],brr[ir]);
          }

          if(arr[m+1]>arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
            SWAP(brr[m+1],brr[ir]);
          }

          if(arr[m]>arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
            SWAP(brr[m],brr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          b=brr[m+1];
          for(;;)
          {
             do i++; while (arr[i]<a);
             do j--; while (arr[j]>a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
             SWAP(brr[i],brr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
          brr[m+1]=brr[j];
          brr[j]=b;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}

/*   function to sort an array into a decreasing order:  a>b>c>....    */  
/* =============================== */
void ArraySort2_Int2(int n, int *arr, int *brr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,NSTACK=50,istack[NSTACK];
     int a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             b=brr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]>=a) break;
                arr[i+1]=arr[i];
                brr[i+1]=brr[i];
             }
             arr[i+1]=a;
             brr[i+1]=b;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);
          SWAP(brr[k],brr[m+1]);

          if(arr[m]<arr[ir])
          {
            SWAP(arr[m],arr[ir]);
            SWAP(brr[m],brr[ir]);
          }

          if(arr[m+1]<arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
            SWAP(brr[m+1],brr[ir]);
          }

          if(arr[m]<arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
            SWAP(brr[m],brr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          b=brr[m+1];
          for(;;)
          {
             do i++; while (arr[i]>a);
             do j--; while (arr[j]<a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
             SWAP(brr[i],brr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
          brr[m+1]=brr[j];
          brr[j]=b;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}

/* =============================== */
void ArraySort_Mix3(int n, B64_long *arr, int *brr, int *crr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,c,NSTACK=50,istack[NSTACK];
     B64_long a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             b=brr[j];
             c=crr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]<=a) break;
                arr[i+1]=arr[i];
                brr[i+1]=brr[i];
                crr[i+1]=crr[i];
             }
             arr[i+1]=a;
             brr[i+1]=b;
             crr[i+1]=c;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);
          SWAP(brr[k],brr[m+1]);
          SWAP(crr[k],crr[m+1]);

          if(arr[m]>arr[ir])
          {
            SWAP(arr[m],arr[ir]);
            SWAP(brr[m],brr[ir]);
            SWAP(crr[m],crr[ir]);
          }

          if(arr[m+1]>arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
            SWAP(brr[m+1],brr[ir]);
            SWAP(crr[m+1],crr[ir]);
          }

          if(arr[m]>arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
            SWAP(brr[m],brr[m+1]);
            SWAP(crr[m],crr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          b=brr[m+1];
          c=crr[m+1];
          for(;;)
          {
             do i++; while (arr[i]<a);
             do j--; while (arr[j]>a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
             SWAP(brr[i],brr[j]);
             SWAP(crr[i],crr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
          brr[m+1]=brr[j];
          brr[j]=b;
          crr[m+1]=crr[j];
          crr[j]=c;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}


/*   to swap the string arrays           */
/* ============================================= */
void s_swap(char **Pair_Name, int i, int j)
/* ============================================= */
{
     char temp[Max_N_NameBase];

     strcpy(temp,Pair_Name[i]);
     strcpy(Pair_Name[i],Pair_Name[j]);
     strcpy(Pair_Name[j],temp);
}


/*   to sort the string array in order          */
/* ============================================= */
void ArraySort_String(int n, char **Pair_Name, int *brr)
/* ============================================= */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,NSTACK=50,istack[NSTACK];
     int temp,MIN=7;
     char p[Max_N_NameBase];

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             strcpy(p,Pair_Name[j]);
             b=brr[j];
             for(i=j-1;i>=m;i--)
             {
                if(strcmp(Pair_Name[i],p)<=0) break;
                strcpy(Pair_Name[i+1],Pair_Name[i]);
                brr[i+1]=brr[i];
             }
             strcpy(Pair_Name[i+1],p);
             brr[i+1]=b;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          s_swap(Pair_Name,k,m+1);
          SWAP(brr[k],brr[m+1]);

          if(strcmp(Pair_Name[m],Pair_Name[ir])>0)
          {
            s_swap(Pair_Name,m,ir);
            SWAP(brr[m],brr[ir]);
          }

          if(strcmp(Pair_Name[m+1],Pair_Name[ir])>0)
          {
            s_swap(Pair_Name,m+1,ir);
            SWAP(brr[m+1],brr[ir]);
          }

          if(strcmp(Pair_Name[m],Pair_Name[m+1])>0)
          {
            s_swap(Pair_Name,m,m+1);
            SWAP(brr[m],brr[m+1]);
          }

          i=m+1;
          j=ir;
          strcpy(p,Pair_Name[m+1]);
          b=brr[m+1];
          for(;;)
          {
             do i++; while (strcmp(Pair_Name[i],p)<0);
             do j--; while (strcmp(Pair_Name[j],p)>0);
             if(j<i) break;
             s_swap(Pair_Name,i,j);
             SWAP(brr[i],brr[j]);
          }
          strcpy(Pair_Name[m+1],Pair_Name[j]);
          strcpy(Pair_Name[j],p);
          brr[m+1]=brr[j];
          brr[j]=b;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}

/*   to swap the string arrays           */
/* ============================================= */
void s_swap2(char **Pair_Name, int i, int j)
/* ============================================= */
{
     char temp[Max_N_NameBase2];

     strcpy(temp,Pair_Name[i]);
     strcpy(Pair_Name[i],Pair_Name[j]);
     strcpy(Pair_Name[j],temp);
}


/*   to sort the string array in order          */
/* ============================================= */
void ArraySort_String2(int n, char **Pair_Name, int *brr)
/* ============================================= */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,NSTACK=50,istack[NSTACK];
     int temp,MIN=7;
     char p[Max_N_NameBase2];

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             strcpy(p,Pair_Name[j]);
             b=brr[j];
             for(i=j-1;i>=m;i--)
             {
                if(strcmp(Pair_Name[i],p)<=0) break;
                strcpy(Pair_Name[i+1],Pair_Name[i]);
                brr[i+1]=brr[i];
             }
             strcpy(Pair_Name[i+1],p);
             brr[i+1]=b;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          s_swap2(Pair_Name,k,m+1);
          SWAP(brr[k],brr[m+1]);

          if(strcmp(Pair_Name[m],Pair_Name[ir])>0)
          {
            s_swap2(Pair_Name,m,ir);
            SWAP(brr[m],brr[ir]);
          }

          if(strcmp(Pair_Name[m+1],Pair_Name[ir])>0)
          {
            s_swap2(Pair_Name,m+1,ir);
            SWAP(brr[m+1],brr[ir]);
          }

          if(strcmp(Pair_Name[m],Pair_Name[m+1])>0)
          {
            s_swap2(Pair_Name,m,m+1);
            SWAP(brr[m],brr[m+1]);
          }

          i=m+1;
          j=ir;
          strcpy(p,Pair_Name[m+1]);
          b=brr[m+1];
          for(;;)
          {
             do i++; while (strcmp(Pair_Name[i],p)<0);
             do j--; while (strcmp(Pair_Name[j],p)>0);
             if(j<i) break;
             s_swap2(Pair_Name,i,j);
             SWAP(brr[i],brr[j]);
          }
          strcpy(Pair_Name[m+1],Pair_Name[j]);
          strcpy(Pair_Name[j],p);
          brr[m+1]=brr[j];
          brr[j]=b;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}


/* creat an int matrix with subscript ange m[nrl...nrh][ncl...nch]  */
int     **imatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch)
{
        B64_long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
        int  **m;

        /* allocate pointers to rows        */
        if((m=(int **)calloc(nrow,sizeof(int*)))==NULL)
        {
           printf("error imatrix: calloc error No. 1 \n");
           return(NULL);
        }
        m+=0;
        m-=nrl;

        /* allocate rows and set pointers to them        */
        if((m[nrl]=(int *)calloc(nrow*ncol,sizeof(int)))==NULL)
        {
           printf("error imatrix: calloc error No. 2 \n");
           return(NULL);
        }
        m[nrl]+=0;
        m[nrl]-=nrl;

        for(i=nrl+1;i<=nrh;i++)
           m[i]=m[i-1]+ncol;
        /* return pointer to array of pointers to rows   */
        return m;
}

/* creat char matrix with subscript ange cm[nrl...nrh][ncl...nch]  */
char    **cmatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch)
{
        B64_long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
        char **cm;

        /* allocate pointers to rows        */
        if((cm=(char **)calloc(nrow,sizeof(char*)))==NULL)
        {
           printf("error cmatrix: calloc error No. 1 \n");
           return(NULL);
        }
        cm+=0;
        cm-=nrl;

        /* allocate rows and set pointers to them        */
        if((cm[nrl]=(char *)calloc(nrow*ncol,sizeof(char)))==NULL)
        {
           printf("error cmatrix: calloc error No. 2 \n");
           return(NULL);
        }
        cm[nrl]+=0;
        cm[nrl]-=nrl;

        for(i=nrl+1;i<=nrh;i++)
           cm[i]=cm[i-1]+ncol;
        /* return pointer to array of pointers to rows   */
        return cm;
}

