import numpy as np
from matplotlib import pyplot
import pandas as pd

#Reading from file

filename = "example.txt" #Change the filename here

pca_file = np.genfromtxt(filename,dtype=str,delimiter="\t")

pca_file = np.matrix(pca_file)
(rows,cols) = pca_file.shape

#Determining the number of features
a = [x for x in range(cols-1)]

required_cols = tuple(a)

pca_features = np.loadtxt(filename, usecols=required_cols)
pca_class = np.loadtxt(filename,dtype=str, usecols=(cols-1,))

mean = np.mean(pca_features,axis=0)
adjusted_matrix = pca_features - mean


#Covariance of adjusted matrix
covariance = np.cov(adjusted_matrix,rowvar=False)

#print covariance.shape (4,4)

#Finding eigen values and eigen vectors
eigen_val , eigen_vector = np.linalg.eig(covariance)

#sorted in ascending order
args = eigen_val.argsort()
#print args

#Determining the indices of the new dimensions.
dim1 = cols-2
dim2 = cols-3

dim1_index= args[dim1]
dim2_index = args[dim2]


eigen_vector1 = eigen_vector[dim1_index]
eigen_vector2 = eigen_vector[dim2_index]

eigen_vector1 = np.matrix(eigen_vector1)
eigen_vector1 = eigen_vector1.transpose()

eigen_vector2 = np.matrix(eigen_vector2)
eigen_vector2 = eigen_vector2.transpose()

#Determing the values of the new dimensions.
new_dimension1 = adjusted_matrix * eigen_vector1
new_dimension2 = adjusted_matrix * eigen_vector2

new_dimension1 = new_dimension1.tolist()
new_dimension2 = new_dimension2.tolist()

new_dimension_data = np.column_stack((new_dimension1,new_dimension2,pca_class))

print "New data"
print new_dimension_data

#Getting the unique classes to determine the colors.
classes = set(pca_class)
classes_list = list(classes)

print "Number of unique classes: ",
print len(classes_list)

colors = ("red","gold","blue","green","black","yellow","orange","maroon","aqua","plum","seagreen")

#Assigning colors to the classes(disease)
index = 0
for i in classes_list:
    pca_class[pca_class == i] = colors[index]
    index = index+1


new_dimension_data = np.column_stack((new_dimension1,new_dimension2,pca_class))

df = pd.DataFrame(new_dimension_data,columns = ["feature1","feature2","disease"])

df["feature1"] = df["feature1"].astype('float')
df["feature2"] = df["feature2"].astype('float')

#Scatter plot with the new dimensions.
pyplot.title("PCA on "+filename)
for i in range(rows):
    pyplot.scatter(df.feature1, df.feature2, c=df.disease)
pyplot.show(block =True)
