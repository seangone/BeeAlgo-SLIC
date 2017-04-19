# BeeAlgo-SLIC

## Introduction

This project is for my graduation project. It is one application of **Bee Algorithm** in **SLIC** Superpixels Segmentation.


##Abstract

First of all, the basic bee algorithm is summarized. The biological basis and model of bee algorithm are introduced. Then, the basic advantages and disadvantages of the algorithm is analyzed. Secondly, the superpixel segmentation problem and the application of the SLIC method are introduced in this paper. Finally, bee algorithm is experimentally applied to the SLIC superpixels algorithm. Compared with other similar application, the innovation of the paper is that the superpixel density of certain amount of pixels of the area instead of the centers of superpixels is taken directly as the parameters of the location of nectar source, so the swarm algorithm does not need to search the center of pixels directly, but is based on the automatic update of the center of the SLIC method. It mainly solves the problems of computing resource scheduling. For those applications who have constraint in the number of super pixels and relatively less constraint in regularity of the superpixels, the algorithm can well realize segmentation requirement: with the same amount of superpixels achieved, it achieve better edge joint. Or with same computing time, there will be a great increase in computing performance. In addition, the method can be applied to application with limited computing resource and the need to focus on some specific areas. However the algorithm is fairly immature, the deficiencies and among them that can be improved are also analyzed.

## Keywords
**bee algorithm**, **superpixels**, **segmentation**, **SLIC**, **resource allocation**

## Files
Here, I give the source of the code. If you are interested in my paper, plase contact me via email Shiyuan AT usc.edu