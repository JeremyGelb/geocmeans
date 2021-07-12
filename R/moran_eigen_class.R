#a concept for using moran spatial vectors in fuzzy clustering.


# In the SFCM, we combine two matrices, the semantic one and a spatial
# one which is just the lagged version of the first.
# I wonder if I could replace it with the spatial eigenvectors. To do so
# I need first to remove all the non relevant vectors (with a thresold value on moran I)
# then I need to find a way to reweight everything because a lot of vectors
# will have a too eavy impact on classification quality
# to do so, I propose the following :
#    1. calculating the overall inertia of the semantic matrix : oI (the sum of the distance matrix)
#    2. calculating the distance matrix using the selected eigen vectors
#    3. rescale that matrix so that it sum is Oi
#    4. re-calculate the coordinates of the vectors with the formula here : https://stackoverflow.com/questions/18096783/using-distance-matrix-to-find-coordinate-points-of-set-of-points
# and here : https://stackoverflow.com/questions/10963054/finding-the-coordinates-of-points-from-distance-matrix
#    5. use this column in the calculus for the spatial matrix, alpha is still to be defined by the user.
# the idea is thus to integrate the spatial eigen vectors but to limit their contribution to a given number of times the ones of the dataset.


# for raster support :
# https://rdrr.io/bioc/OLIN/man/ma.matrix.html
# https://www.rdocumentation.org/packages/raster/versions/3.4-10/topics/focal
