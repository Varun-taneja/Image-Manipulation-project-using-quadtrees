#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "proj.h"




void freeMatrix(colour_p **mat, int H){   // releases memeory ie frees memory
	int i;
	for(i = 0; i < H; i++){
		free(mat[i]);
	}
	free(mat);
}

void freeQT(node *root){     // Memory release for a tree
	if(root == NULL)
		return;

	freeQT(root->NW);
	freeQT(root->NE);
	freeQT(root->SE);
	freeQT(root->SW);

	free(root);
	root = NULL;
}

void nodeAdd(node **root, colour_p p, int l){  // addition of a node in a quaternary tree
	node *q = (node*)malloc(sizeof(node));
	q->R = p.R;
	q->G = p.G;
	q->B = p.B;
	q->A = l*l;
	q->NW = NULL;
	q->NE = NULL;
	q->SE = NULL;
	q->SW = NULL;

	if(*root == NULL){
		*root = q;
		return;
	}
	else{
		if((*root)->NW == NULL){
			(*root)->NW = q;
		}
		else if((*root)->NE == NULL){
			(*root)->NE = q;
		}
		else if((*root)->SE == NULL){
			(*root)->SE = q;
		}
		else if((*root)->SW == NULL){
			(*root)->SW = q;
		}
	}
}

int Uniformity(colour_p **mat, int x, int y, int l, colour_p *p){    // calculating the degree of uniformity of in area of ​​the matrix and the average color
	long long avg = 0;
	int i, j;
	long long R, G, B;
	long long r, g, b;
	long long size;
	size = l;
	R = 0;
	G = 0;
	B = 0; 	
	r = g = b =0;

	for(i = x; i < x + l; i++)
		for(j = y; j < y + l; j++){	
			r= mat[i][j].R;
			R = R + r;
			
			g = mat[i][j].G;
			G = G + g;
			
			b = mat[i][j].B;
			B = B + b;	
		}
		
		R = R/(size * size);
		G = G/(size * size);
		B = B/(size * size);	


		(*p).R = R;
		(*p).G = G;
		(*p).B = B;
	
	for(i = x; i < x + l; i++)
		for(j = y; j < y + l; j++)
			{
				avg += (R - mat[i][j].R)*(R - mat[i][j].R);
				avg += (B - mat[i][j].B)*(B - mat[i][j].B);
				avg += (G - mat[i][j].G)*(G - mat[i][j].G);
			}
	
	avg = avg/(3 * l * l);
	return avg;
}

void Compression(node **root, int x, int y, int l, int W, colour_p **mat, int bound, int *c_num, int *n_num){	// Compression in a quaternary tree of pixels that have RGB values 
	int average;  //average - current level of uniformity   
	colour_p p;   //p - the average color for the current square
	average = Uniformity(mat, x, y, l, &p);

	if(l < 1)
		return;

	(*n_num) = (*n_num) + 1;
	nodeAdd(root, p, l);

	if(average > bound){
		Compression(&(*root)->NW, x, y, l/2, W, mat,bound, c_num, n_num);
		(*c_num) = (*c_num)+1;
		Compression(&(*root)->NE, x, y + l/2, l/2, W, mat,bound, c_num, n_num);
		(*c_num) = (*c_num)+1;
		Compression(&(*root)->SE, x + l/2, y + l/2, l/2, W, mat,bound, c_num, n_num);
		(*c_num) = (*c_num)+1;
		Compression(&(*root)->SW, x + l/2, y, l/2, W, mat, bound, c_num, n_num);	
		(*c_num) = (*c_num)+1;
	}	
}

void makeImage(colour_p **mat, int W, int H, int max_num, char *file_name){    // making an image from an array of pixels that have RGB values
	
	FILE *f = fopen(file_name,"wb");

	fprintf(f, "%s\n", "P6");
	fprintf(f, "%d %d\n", H, W);
    fprintf(f, "%d\n", max_num);
	
    int i,j;
  
	for(i = 0; i < H; i++){
	 	fwrite(mat[i],sizeof(colour_p),W,f);
	}
	fclose(f);
}

int powTwo(int x){       //function for obtaining a number ie a power of 2
	int two = 2;
	while(two * two != x)
		two = two * 2;
	return two;
}

int Depth(node* root){    // determining the depth of a quaternary tree
    if (root == NULL)
        return 0;
    else{       
        int lf = Depth(root->NW);
        int lr = Depth(root->NE);
        int br = Depth(root->SE);
        int bl = Depth(root->SW); 
       
    	int max = lf;
    	if(lr > max)
    		max = lr;
    	if(br > max)
    		max = br;
    	if(bl > max)
    		max = bl;
    	return(max+1);
    }
}

void Level(QTNode *v, int *k, int *child_num, node *root, int level){     // function for loading a QuadNode array by traversing a tree in levels
    if (root == NULL)
        return;
   
    if (level == 1){	
        v[*k].R = root->R;
        v[*k].G = root->G;
        v[*k].B = root->B;
        v[*k].A = root->A;

        if(root->NW == NULL){
        	v[*k].NW = -1;
        	v[*k].NE = -1;
        	v[*k].SE = -1;
        	v[*k].SW = -1;
        }
        else{	
        	v[*k].NW = *child_num;
        	*child_num = (*child_num)+1;
        	v[*k].NE = *child_num;
        	*child_num = (*child_num)+1;
        	v[*k].SE = *child_num;
        	*child_num = (*child_num)+1;
        	v[*k].SW = *child_num;
        	*child_num = (*child_num)+1;
        }
        *k = (*k)+1;
    }
    else if (level > 1){
        Level(v, k, child_num, root->NW, level-1);
        Level(v, k, child_num, root->NE, level-1);
        Level(v, k, child_num, root->SE, level-1);
        Level(v, k, child_num, root->SW, level-1);
    }
}

// bringing a quaternary tree to a shape in which each level contains a fixed 4 ^ level ,the root is consider at level 0 in this case
void Complete_QT(node **root, int l, int depth, int line) {	
	if(*root == NULL)
		return;

	if(line < depth && (*root)->NW == NULL){	
		colour_p p;
		p.R = (*root)->R;
		p.G = (*root)->G;
		p.B = (*root)->B;

		nodeAdd(root,p,l);
		nodeAdd(root,p,l);
		nodeAdd(root,p,l);
		nodeAdd(root,p,l);
	}

	Complete_QT(&(*root)->NW,l/2,depth,line+1);
	Complete_QT(&(*root)->NE,l/2,depth,line+1);
	Complete_QT(&(*root)->SE,l/2,depth,line+1);
	Complete_QT(&(*root)->SW,l/2,depth,line+1);
}

void fileOut(QTNode *v, int *k, int *child_num, node *root, int depth, char *file_name, int c_num, int n_num, int W){    //creating a Comprressed file from a given array
	int i;	

	for(i = 1; i <= depth; i++){
		Level(v, k, child_num, root, i);
	}
	FILE *f = fopen(file_name, "wb");
	
	fwrite(&c_num, sizeof(int), 1, f);
	fwrite(&n_num, sizeof(int), 1, f);	
	
	for(i = 0; i < *k; i++){	
		fwrite(&v[i], sizeof(QTNode), 1, f);
	}		
	fclose(f);
}

void readBinary(QTNode *v, int *c_num, int *n_num, char *file_name){    // Reading an array from the compressed file
	FILE *f = fopen(file_name, "rb");

	fread(c_num, sizeof(int), 1, f);
	fread(n_num, sizeof(int), 1, f);

	int i;

	for(i = 0; i < *n_num; i++){
		fread(&v[i].R, sizeof(unsigned char), 1, f);
		fread(&v[i].G, sizeof(unsigned char), 1, f);
		fread(&v[i].B, sizeof(unsigned char), 1, f);
		fread(&v[i].A, sizeof(int), 1, f);
		fread(&v[i].NW, sizeof(int), 1, f);
		fread(&v[i].NE, sizeof(int), 1, f);
		fread(&v[i].SW, sizeof(int), 1, f);
		fread(&v[i].SE, sizeof(int), 1, f);		
	}
}

void matValues(colour_p **mat, int x, int y, int l,unsigned char R, unsigned char G, unsigned char B){   // load RGB values on a pixel in a square of the matrix
	int i,j;
	for(i = x; i < x+l; i++)
		for(j = y; j < y+l; j++){
				mat[i][j].R = R;
				mat[i][j].G = G;
				mat[i][j].B = B;
		}
}

void loadQT(node **root, QTNode *v, int l, int i, int n_num){    // Load a vector into a quaternary tree
	if(n_num == 0 || l < 1)
		return;

	colour_p p;
	p.R = v[i].R;
	p.G = v[i].G;
	p.B = v[i].B;

	nodeAdd(root,p,l);

	if(v[i].NW + v[i].NE + v[i].SE + v[i].SW != -4 ){
		loadQT(&(*root)->NW, v, l/2, v[i].NW, n_num - 1);
		loadQT(&(*root)->NE, v, l/2, v[i].NE, n_num - 1);
		loadQT(&(*root)->SE, v, l/2, v[i].SE, n_num - 1);
		loadQT(&(*root)->SW, v, l/2, v[i].SW, n_num - 1);
	}	
}

void loadMat(node *root, colour_p **mat, int x, int y, int l, int *c_num, int n_num){	// Loading a matrix with values ​​from a  quadtree
		
	if((*c_num) == 0 || l == 0 )
		return;

	if(root->NW == NULL){	 
	 	matValues(mat, x, y, l, root->R, root->G, root->B); 
	 	*c_num = (*c_num) - 1;
	}
	else{		 	
	 	loadMat(root->NW, mat, x, y, l/2, c_num, n_num);
	 	loadMat(root->NE, mat, x, y + l/2, l/2, c_num, n_num);
	 	loadMat(root->SW, mat, x + l/2, y, l/2, c_num, n_num);
		loadMat(root->SE, mat, x + l/2, y + l/2, l/2, c_num, n_num);	 	
	}
}

//Loading a matrix with the "averages" of colors on the leaves of two quaternary trees with the same number of levels and nodes
void Overlap(node *root1, node *root2, colour_p **mat, int x, int y, int l){		
	if( l == 0 )
		return;

	if(root1->NW == NULL){	 
	 	matValues(mat, x, y, l, (root1->R+root2->R)/2, (root1->G + root2->G)/2, (root1->B + root2->B)/2); 
	}
	else{	 	
	 	Overlap(root1->NW, root2->NW, mat, x, y, l/2);
	 	Overlap(root1->NE, root2->NE, mat, x, y + l/2, l/2);
	 	Overlap(root1->SW, root2->SW, mat, x + l/2, y, l/2);
		Overlap(root1->SE, root2->SE, mat, x + l/2, y + l/2, l/2);
	}
}

// making the whole image of only red colour
void red(node *root1, colour_p **mat, int x, int y, int l){		
	if( l == 0 )
		return;

	if(root1->NW == NULL){	 
	 	matValues(mat, x, y, l, root1->R, 0, 0); 
	}
	else{	 	
	 	red(root1->NW, mat, x, y, l/2);
	 	red(root1->NE, mat, x, y + l/2, l/2);
	 	red(root1->SW, mat, x + l/2, y, l/2);
		red(root1->SE, mat, x + l/2, y + l/2, l/2);
	}
}

// making the whole image of only green colour
void green(node *root1, colour_p **mat, int x, int y, int l){		
	if( l == 0 )
		return;

	if(root1->NW == NULL){	 
	 	matValues(mat, x, y, l, 0, root1->G, 0); 
	}
	else{	 	
	 	green(root1->NW, mat, x, y, l/2);
	 	green(root1->NE, mat, x, y + l/2, l/2);
	 	green(root1->SW, mat, x + l/2, y, l/2);
		green(root1->SE, mat, x + l/2, y + l/2, l/2);
	}
}

// making the whole image of only blue colour
void blue(node *root1, colour_p **mat, int x, int y, int l){		
	if( l == 0 )
		return;

	if(root1->NW == NULL){	 
	 	matValues(mat, x, y, l, 0, 0, root1->B); 
	}
	else{	 	
	 	blue(root1->NW, mat, x, y, l/2);
	 	blue(root1->NE, mat, x, y + l/2, l/2);
	 	blue(root1->SW, mat, x + l/2, y, l/2);
		blue(root1->SE, mat, x + l/2, y + l/2, l/2);
	}
}
/* making the whole image into greysacle by setting the RGB values for a pixel to 0.3*r + 0.59*g + 0.11*b */
void greyscale(node *root1, colour_p **mat, int x, int y, int l){		
	if( l == 0 )
	return;

 	if(root1->NW == NULL){	 
 		int c = (root1->R * 0.3) + (root1->G * 0.59) + (root1->B * 0.11);
 		matValues(mat, x, y, l, c, c, c);  		
 	}
 	else{	 	
 	 	greyscale(root1->NW, mat, x, y, l/2);
 	 	greyscale(root1->NE, mat, x, y + l/2, l/2);
 	 	greyscale(root1->SW, mat, x + l/2, y, l/2);
 		greyscale(root1->SE, mat, x + l/2, y + l/2, l/2);
 	}
}

//Making the whole picture black and white (255,255,255) = white , (0,0,0) = black
void bw(node *root1, colour_p **mat, int x, int y, int l){		
	if( l == 0 )
		return;

	if(root1->NW == NULL){	 
		int c = (root1->R + root1->G + root1->B)/3 ;
		if(c > 200){
	 		matValues(mat, x, y, l, 0, 0, 0); 
		}
		else
			matValues(mat, x, y, l, 255, 255, 255); 	
	}
	else{	 	
	 	bw(root1->NW, mat, x, y, l/2);
	 	bw(root1->NE, mat, x, y + l/2, l/2);
	 	bw(root1->SW, mat, x + l/2, y, l/2);
		bw(root1->SE, mat, x + l/2, y + l/2, l/2);
	}
}

// Mirroring an image horizontally
void Horizontal(node **root){
	if(*root == NULL)
		return;

	if((*root)->NW != NULL){
		node *temp = NULL;
		
		temp = (*root)->NW; //swap NW with NE
		(*root)->NW = (*root)->NE;
		(*root)->NE = temp;
		
		temp = (*root)->SW; //swap SW with SE
		(*root)->SW = (*root)->SE;
		(*root)->SE = temp;

		Horizontal(&(*root)->NW);
		Horizontal(&(*root)->NE);
		Horizontal(&(*root)->SE);
		Horizontal(&(*root)->SW);
	}
}
// Mirroring an image vertically 
void Vertical(node **root){
	if(*root == NULL)
		return;

	if((*root)->NW != NULL){	
		node *temp = NULL;

		temp = (*root)->NW; //swap NW with SW
		(*root)->NW = (*root)->SW;
		(*root)->SW = temp;
		
		temp = (*root)->NE; //swap NE with SE
		(*root)->NE = (*root)->SE;
		(*root)->SE = temp;
		
		Vertical(&(*root)->NW);
		Vertical(&(*root)->NE);
		Vertical(&(*root)->SE);
		Vertical(&(*root)->SW);	
	}
}

