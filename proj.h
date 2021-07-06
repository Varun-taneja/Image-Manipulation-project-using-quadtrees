typedef struct QTNode {
	unsigned char R,G,B;
	unsigned int A;
	int NW, NE;
	int SW, SE;
}__attribute__((packed)) QTNode;

typedef struct colour_p{
	unsigned char R;
	unsigned char G;
	unsigned char B;
}colour_p;

typedef struct node{
	unsigned char R, G, B;
	unsigned int A;
	struct node *NW;
	struct node *NE;
	struct node *SE;
	struct node *SW;
}node;

void freeMatrix(colour_p **mat, int H);
void freeQT(node *root);
void nodeAdd(node **root, colour_p p, int l);
int Uniformity(colour_p **mat, int x, int y, int l, colour_p *p);
void Compression(node **root, int x, int y, int l, int W, colour_p **mat, int bound, int *c_num, int *n_num);
void makeImage(colour_p **mat, int W, int H, int max_num, char *file_name);
int powTwo(int x);
int Depth(node* root);
void Level(QTNode *v, int *k,int *child_num, node *root, int level);
void Complete_QT(node **root, int l, int depth, int line);
void fileOut(QTNode *v, int *k, int *child_num, node *root, int depth, char *file_name, int c_num, int n_num, int W);
void readBinary(QTNode *v, int *c_num, int *n_num, char *file_name);
void matValues(colour_p **mat, int x, int y, int l,unsigned char R, unsigned char G, unsigned char B);
void loadQT(node **root, QTNode *v, int l, int i, int n_num);
void loadMat(node *root, colour_p **mat, int x, int y, int l, int *c_num, int n_num);
void Overlap(node *root1, node *root2, colour_p **mat, int x, int y, int l);
void red(node *root1, colour_p **mat, int x, int y, int l);
void green(node *root1, colour_p **mat, int x, int y, int l);
void blue(node *root1, colour_p **mat, int x, int y, int l);
void greyscale(node *root1, colour_p **mat, int x, int y, int l);
void bw(node *root1, colour_p **mat, int x, int y, int l);
void Horizontal(node **root);
void Vertical(node **root);

