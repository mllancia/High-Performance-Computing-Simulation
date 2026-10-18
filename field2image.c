/* field2image accepts a sequence of configurations and makes
 * a images or a movie
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <err.h>
#include <math.h>
#include <getopt.h>
#include <dirent.h>
#include <png.h>

#ifndef M_PI
#define M_PI 3.14159265358979323844
#endif

enum {
    max_file_list=65536,
    max_file_name=512,
};

struct entry {
    char name[max_file_name];
    int seq;
};

struct entry File_list[max_file_list];
int File_n;
int Lx,Ly; /* lattice size */

void read_directory(const char *name)/*{{{*/
{
    DIR *D=opendir(name);
    struct dirent *E;
    if(D==NULL) err(1,"Could not read the directory %s",name);
    while((E=readdir(D))) {
        int seq=0;
        char extra;
        if(sscanf(E->d_name,"seq.%d%c",&seq,&extra)==1) {
            if(File_n>=max_file_list) errx(1,"Too many input files\n");
            snprintf(File_list[File_n].name,
                    sizeof(File_list[File_n].name),"%s/%s",name,E->d_name);
            File_list[File_n].seq=seq;
            File_n++;
        }
    }
    closedir(D);
}/*}}}*/

int cmp_entry(const void *a,const void *b)/*{{{*/
{
    const struct entry *A=a;
    const struct entry *B=b;
    if(A->seq==B->seq) errx(1,"Two files with the same sequence number %s and %s\n",A->name,B->name);
    return A->seq-B->seq;
}/*}}}*/

void sort_list()/*{{{*/
{
    qsort(File_list,File_n,sizeof(File_list[0]),cmp_entry);
}/*}}}*/

int make_png_image(const char *file,const float *red,const float *green,const float *blue,int nx,int ny)/*{{{*/
{
    FILE *F=fopen(file,"w");
    int error=0;
    if(!F) err(1,"Could not create %s",file);
    png_structp P=png_create_write_struct(PNG_LIBPNG_VER_STRING,
            NULL,NULL,NULL);
    if(P==NULL) { error=1; goto close;}
    png_infop  I=png_create_info_struct(P);
    if(I==NULL) { error=2; goto png_destroy; }
    if (setjmp(png_jmpbuf(P))) {
        png_destroy_write_struct(&P, &I);
        goto close;
    }
   png_init_io(P,F);
   png_set_IHDR(P,I,nx,ny,8,PNG_COLOR_TYPE_RGB,PNG_INTERLACE_NONE,PNG_COMPRESSION_TYPE_DEFAULT,PNG_FILTER_TYPE_DEFAULT);
   png_write_info(P,I);
   {
       png_byte *row_pointers[ny];
       unsigned char *image_array=malloc(nx*ny*3);
       if(!image_array) errx(1,"Could not allocated RGB char arrays\n");
       for(int y=0;y<ny;y++) 
           row_pointers[y]=&image_array[y*nx*3];
       for(int y=0,i=0;y<ny;y++) 
           for(int x=0;x<nx;x++,i++) {
               int r=255*red[i];   if(r<0) r=0; if(r>255) r=255;
               int g=255*green[i]; if(g<0) g=0; if(g>255) g=255;
               int b=255*blue[i];  if(b<0) b=0; if(b>255) b=255;
               image_array[y*nx*3+3*x+0]=r;
               image_array[y*nx*3+3*x+1]=g;
               image_array[y*nx*3+3*x+2]=b;
           }
       png_write_image(P, row_pointers);
       free(image_array);
   }
   png_write_end(P,I);
   png_destroy_write_struct(&P,&I);
png_destroy:
   png_destroy_write_struct(&P,(png_infopp)NULL);
close:
   fclose(F);
   return error;
}/*}}}*/

void make_image(struct entry *E)/*{{{*/
{
    FILE *F=fopen(E->name,"r");
    if(F==NULL) errx(1,"Could not read %s",E->name);
    int latsize[4];
    float *data;
    if(fread(latsize,sizeof(latsize[0]),4,F)!=4 
            || latsize[0]!=826886224
            || latsize[1]!=1
            || latsize[2]<2 || latsize[3]<2 
            || latsize[2]>65536 || latsize[3]>65536 
      ) goto bad_format;
    if(Ly && latsize[2]!=Ly) goto incompatible;
    if(Lx && latsize[3]!=Lx) goto incompatible;
    Ly=latsize[2];
    Lx=latsize[3];
    data=malloc(sizeof(data[0])*latsize[0]*latsize[1]);
    if(data==NULL) goto allocation;
    if(fread(data,sizeof(data[0]),Lx*Ly,F)!=Lx*Ly) goto bad_format;
    char pngfile[256];
    if(snprintf(pngfile,sizeof(pngfile),"%s.png",E->name)>=sizeof(pngfile))
        goto long_png_name;
    int ls=Lx*Ly;
    {
        float *red=malloc(sizeof(float)*ls);
        float *green=malloc(sizeof(float)*ls);
        float *blue=malloc(sizeof(float)*ls);
        if(!red || !green || !blue) errx(1,"RGB array allocation failure");
//#pragma omp parallel for
        for(int i=0;i<ls;i++) {
            double phi=tanh(data[i])*0.5+0.5;
            red[i]=sin(0.5*M_PI*phi);
            green[i]=pow(2*fabs(phi-0.5),2);
            blue[i]=cos(0.5*M_PI*phi);
        }
        make_png_image(pngfile,red,green,blue,Lx,Ly);
        free(red);
        free(green);
        free(blue);
    }
    free(data);
    fclose(F);
    return;
bad_format:
    errx(1,"File %s was in bad format\n",E->name);
allocation:
    errx(1,"Memory allocation error while processing %s\n",E->name);
incompatible:
    errx(1,"This file is incompatible with the others %s\n",E->name);
long_png_name:
    errx(1,"Too long name for the png file derived from %s\n",E->name);
}/*}}}*/

void process_dir(const char *name)/*{{{*/
{
    char cmd[1024];
    read_directory(name);
    sort_list();
#pragma omp parallel for schedule(dynamic)
    for(int i=0;i<File_n;i++) {
        puts(File_list[i].name);
        make_image(&File_list[i]);
    }
    if(sizeof(cmd)-1<
    snprintf(cmd,sizeof(cmd),
            "ffmpeg -pattern_type glob -i '%s/seq.*.png' %s.flv",name,name))
        errx(1,"Too long command name for ffmpeg\n");
    if (system(cmd) != 0) {
        warnx("system command failed: %s", cmd);
    }
    printf("The result is in %s.flv\n",name);
}/*}}}*/

int main(int argc,char **argv)/*{{{*/
{
    const char *Help="field2image\n"
        "--------------------------------------------------------------\n"
        " field2image directory\n"
        "       will look for seq.0000, seq.0001, seq.0002 ...\n"
        "       files in the given directory \n"
        "       and convert each of them to png images.\n"
        "       Finally, it produces a directory.flv movie file\n";
    int o;
    while((o=getopt(argc,argv,"h"))!=-1) switch(o) {
        case 'h': puts(Help); exit(0);
        default: puts("use -h for help\n"); exit(0);
    }
    for(;optind<argc;optind++) process_dir(argv[optind]);
}/*}}}*/
