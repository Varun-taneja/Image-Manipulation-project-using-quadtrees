# Image-Manipulation-project-using-quadtrees

## What are Quadtrees?
A quadtree is a tree data structure in which each internal node has exactly four children. Quadtrees are the two-dimensional analog of octrees and are most often used to partition a two-dimensional space by recursively subdividing it into four quadrants or regions. The data associated with a leaf cell varies by application, but the leaf cell represents a "unit of interesting spatial information“.


![image](https://user-images.githubusercontent.com/69303551/127506441-b381175f-331c-4222-9afd-1dfc87e59132.png)


## What have I done ?
1)Compression of Images 
2)Decompression of Images
3)Horizontal or Vertical flipping of Images 
4)Overlapping two Images
5)Changing an image into a single colour ie – Red / Green / Blue / greyscale / black&white


![Screenshot (1802)](https://user-images.githubusercontent.com/69303551/127505960-c1cfdd57-606d-49ea-925d-9082ba232953.png)

[NOTE:The Images used in my program are all in PPM (portable pixmap) format.]


## What is a PPM image ?
A PPM file consists of two parts, a header and the image data. The header consists of at least three parts –The first "line" is a magic PPM identifier it can be P6.The next line consists of the width and height of the image as ASCII numbers. The last part of the header gives the maximum value of the colour components for the pixels, this allows the format to describe more than single byte (0 to 255) colour values. 

[NOTE: To view ppm images we have to download FastStone Image Viewer ]

## General Idea -
1)The root of the Quadtree represents a square and the four children of the root correspond to the four squares that making up the whole image and this these 4 squares can then be recursively split into another four children. 
2)Till the time a child node doesn’t have pixels with average colour below the entered boundary, it can’t represent a leaf node and therefore it must be split up into four nodes. 3)We repeat this idea unit all the leaf nodes have pixels with average color that corresponds to the entered boundary.
The average colour of each pixel is calculated using the formula -
![image](https://user-images.githubusercontent.com/69303551/127506046-e27f9cd1-9c68-47ac-8802-0dc03acb9bdf.png)


## How to use -

![Screenshot (1804)](https://user-images.githubusercontent.com/69303551/127506334-744d7841-302a-4f80-adc5-94e6937e1d39.png)
