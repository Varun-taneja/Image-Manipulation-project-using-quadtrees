#include "proj.h" 
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char *argv[]){	
	
	int W, H;
	int max_num;
	int i, j;
	int bound;		
	
	if(strcmp(argv[1], "-c") == 0){          // For compression of an image into an output file
		
		colour_p **mat;
		FILE *f = fopen(argv[3], "rb");
		char type[3];	
		bound = atoi(argv[2]);

		fscanf(f,"%2s",type);
		type[2] ='\0';
		fscanf(f, "%d", &H);
    	fscanf(f, "%d", &W);
    	fscanf(f, "%d", &max_num);
    	char c;
    	fscanf(f, "%c ",&c);
  
    	mat = malloc(H*sizeof(colour_p*));
		for(i = 0; i < H; i++)
			mat[i] = malloc(W*sizeof(colour_p));

		for(i = 0; i < H; i++){
			for(j = 0; j < W; j++){
				fread(&mat[i][j].R, sizeof(unsigned char), 1, f);
			 	fread(&mat[i][j].G, sizeof(unsigned char), 1, f);
			 	fread(&mat[i][j].B, sizeof(unsigned char), 1, f);
			}
		}	
		fclose(f);	

		node *root;	
		root = NULL;
		int c_num = 0, n_num = 0, depth = 0; 
	
		Compression(&root, 0, 0, W, W, mat, bound, &c_num, &n_num);
		depth = Depth(root);

		QTNode *v;
		v = malloc(n_num*sizeof(QTNode));
		int k = 0, child_num = 1;

		fileOut(v, &k, &child_num, root, depth + 1, argv[4], c_num, n_num, W);

		freeMatrix(mat, W);
		if(v != NULL)
			free(v);
		freeQT(root);
	}
	
	
	if(strcmp(argv[1], "-d") == 0){     // decompression of an output file into an image

		colour_p **mat;
		QTNode *v;
		int max_num = 255;
		int c_num, n_num, W;

		FILE *f = fopen(argv[2],"rb");

		fread(&c_num, sizeof(int), 1, f);
		fread(&n_num, sizeof(int), 1, f);
		v = malloc(n_num * sizeof(QTNode));
		fclose(f);
		
		readBinary(v, &c_num, &n_num,argv[2]);

		if(c_num == 0 && n_num == 1)
			c_num = 1;
		
		W = powTwo(v[0].A);	
   	 	mat = malloc(W * sizeof(colour_p*));
		for(i = 0; i < W; i++){
			mat[i] = malloc(W * sizeof(colour_p));
        }

		int k = 0;
		int x = 0, y = 0, l;
		l = W;

		node *tree;	
		tree = NULL;
		
		loadQT(&tree, v, W, k, n_num);

		x = y = 0;
		l = W;

		if(n_num == 1){
			matValues(mat, x, y, l, v[0].R, v[0].G, v[0].B);
		}
		else
			loadMat(tree, mat, x, y, l, &c_num,n_num);
			makeImage(mat, W, W, max_num, argv[3]);
		
		freeMatrix(mat,W);
		if(v != NULL)
			free(v);
		freeQT(tree);
	}	

	
	if(strcmp(argv[1],"-m") == 0){	 // for horizontal and vertical mirrioring
		
		colour_p **mat;
		int mirror;
		if(strcmp(argv[2],"h")==0)
			mirror = 1;  //morror == 1 -> mirror horizonally
		else if(strcmp(argv[2],"v")==0)
			mirror = 2;  //mirror == 2 -> mirror vertically
		else
			mirror = 0;		

		bound = atoi(argv[3]);
		FILE *f = fopen(argv[4],"rb");
		char type[3];
	
		fscanf(f,"%2s",type);
		type[2]='\0';
		fscanf(f, "%d", &H);
    	fscanf(f, "%d", &W);
    	fscanf(f, "%d", &max_num);
    	char c;
    	fscanf(f,"%c",&c);
   
    	mat = malloc(H * sizeof(colour_p*));
		for(i = 0; i < H; i++)
			mat[i] = malloc(W*sizeof(colour_p));

		for(i = 0; i < H; i++){
			for(j = 0; j < W; j++)
				fread(&mat[i][j], sizeof(colour_p), 1, f);
		}
		fclose(f);

		node *root = NULL;
		int c_num = 0, n_num = 0, depth = 0; 
		int x,y,l;

		Compression(&root, 0, 0, W, W, mat, bound, &c_num, &n_num);
		depth = Depth(root);

		if(mirror == 1){
			Horizontal(&root);
		}
		else if (mirror == 2){
			Vertical(&root);
		}
		x = y =0;
		l = W;

		loadMat(root, mat, x, y, l, &c_num,  n_num);
		makeImage(mat, W, W, max_num, argv[5]);
		freeMatrix(mat, W);
		freeQT(root);
	}


	
	if(strcmp(argv[1],"-o") == 0){    // for overlapping 2 images
		
		node *root1 = NULL;
		node *root2 = NULL;
		colour_p **mat1, **mat2;
		int W1, H1, W2, H2;

		bound = atoi(argv[2]);
		FILE *f = fopen(argv[3], "rb");
		char type[3];
	
		fscanf(f,"%2s",type);
		type[2]='\0';
		fscanf(f, "%d", &H1);
		fscanf(f, "%d", &W1);
		fscanf(f, "%d", &max_num);
		char c;
		fscanf(f,"%c",&c);
   
		mat1 = malloc(H1 * sizeof(colour_p*));
		for(i = 0; i < H1; i++)
			mat1[i] = malloc(W1 * sizeof(colour_p));

		for(i = 0; i < H1; i++){
			for(j = 0; j < W1; j++)
				fread(&mat1[i][j], sizeof(colour_p), 1, f);
		}
		fclose(f);

		FILE *g = fopen(argv[4], "rb");		
		
		fscanf(g,"%2s",type);
		type[2]='\0';
		fscanf(g, "%d", &H2);
		fscanf(g, "%d", &W2);
		fscanf(g, "%d", &max_num);
		fscanf(g, "%c", &c);
	
		mat2 = malloc(H2 * sizeof(colour_p*));
		for(i = 0; i < H2; i++)
			mat2[i] = malloc(W2 * sizeof(colour_p));

		for(i = 0; i < H2; i++){
			for(j = 0; j < W2; j++)
				fread(&mat2[i][j], sizeof(colour_p), 1, g);
		}
		fclose(g);

		int c_num1 = 0, n_num1 = 0, depth1 = 0; 
		int c_num2 = 0, n_num2 = 0, depth2 = 0;
		int x = 0, y = 0, l;
		l = W1;

		Compression(&root1, 0, 0, W1, W1, mat1, bound, &c_num1, &n_num1);
		Compression(&root2, 0, 0, W2, W2, mat2, bound, &c_num2, &n_num2);

		depth1 = Depth(root1);
		depth2 = Depth(root2);

		int depth;

		if(depth1 > depth2)
			depth = depth1;
		else
			depth = depth2;		

		int k = 1;

		Complete_QT(&root1, W1, depth, 1);
		Complete_QT(&root2, W2, depth, 1);
		depth1 = Depth(root1);
		depth2 = Depth(root2);
		Overlap(root1, root2, mat1 , x, y, l);
		makeImage(mat1, W1, W1, max_num, argv[5]);

		freeMatrix(mat1, W1);
		freeQT(root1);
		freeMatrix(mat2, W2);
		freeQT(root2);
	}

	if(strcmp(argv[1],"-t") == 0){	 //For changing the original colour of an image
		colour_p **mat;
		int colour;
		if(strcmp(argv[2],"R") == 0)
			colour = 1;  //colour == 1 -> red
		else if(strcmp(argv[2],"G") == 0)
			colour = 2;  //colour == 2 -> green
		else if(strcmp(argv[2],"B")==0)
			colour = 3;  //colour == 3 -> blue
		else if(strcmp(argv[2],"Gr")==0)
			colour = 4;  //colour == 4 -> greyscale
		else if(strcmp(argv[2],"Bw")==0)
			colour = 5;  //colour == 5 -> black and white
		else
			colour = 0;		

		node *root1 = NULL;
		colour_p **mat1;
		int W1, H1;

		bound = atoi(argv[3]);
		FILE *f = fopen(argv[4], "rb");
		char type[3];
	
		fscanf(f, "%2s", type);
		type[2]= '\0';
		fscanf(f, "%d", &H1);
		fscanf(f, "%d", &W1);
		fscanf(f, "%d", &max_num);
		char c;
		fscanf(f,"%c",&c);
   
		mat1 = malloc(H1 * sizeof(colour_p*));
		for(i = 0; i < H1; i++)
			mat1[i] = malloc(W1 * sizeof(colour_p));

		for(i = 0; i < H1; i++){
			for(j = 0; j < W1; j++)
				fread(&mat1[i][j], sizeof(colour_p), 1, f);
		}
		fclose(f);

		int c_num1 = 0, n_num1 = 0, depth1 = 0; 
		int x = 0, y = 0, l;
		l = W1;

		Compression(&root1, 0, 0, W1, W1, mat1, bound, &c_num1, &n_num1);
		depth1 = Depth(root1);		
		Complete_QT(&root1, W1, depth1, 1);

		depth1 = Depth(root1);
		if (colour == 1){
			red(root1, mat1, x, y, l);
		}
		else if (colour == 2){
			green(root1, mat1, x, y, l);
		}
		else if (colour == 3){
			blue(root1, mat1 , x, y, l);
		}
		else if (colour == 4){
			greyscale(root1, mat1 , x, y, l);
		}
		else if (colour == 5){
			bw(root1, mat1 , x, y, l);
		}
		else {
			printf(" wrong input");
			return 0;
		}
		
		makeImage(mat1, W1, W1, max_num, argv[5]);
		freeMatrix(mat1, W1);
		freeQT(root1);
	}
	
	return 0;
}

