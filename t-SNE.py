import numpy as np
from matplotlib import pyplot
import pandas as pd
from sklearn.manifold import TSNE

#Reading from file

filename = "pca_c.txt" #Change the filename here

pca_file = np.genfromtxt(filename,dtype=str,delimiter="\t")

pca_file = np.matrix(pca_file)
(rows,cols) = pca_file.shape

print rows,cols

#Determining the number of features
a = [x for x in range(cols-1)]

required_cols = tuple(a)

pca_features = np.loadtxt(filename, usecols=required_cols)
pca_class = np.loadtxt(filename,dtype=str, usecols=(cols-1,))


#Getting the unique classes to determine the colors.
classes = set(pca_class)
classes_list = list(classes)

print "Number of unique classes: ",
print len(classes_list)

colors = ("red","gold","blue","green","black","yellow","orange","maroon","aqua","plum","seagreen")

#Assigning colors to the classes(disease)
index = 0
for i in classes_list:
    print "index:",index
    pca_class[pca_class == i] = colors[index]
    index = index+1

#Implementing t-SNE on the data.
#Referred from http://scikit-learn.org/stable/modules/generated/sklearn.manifold.TSNE.html
X_embedded = TSNE(n_components=2).fit_transform(pca_features)

print X_embedded.shape

new_dimension1 =  X_embedded[:,0].tolist()
new_dimension2 = X_embedded[:,1].tolist()

new_dimension_data = np.column_stack((new_dimension1,new_dimension2,pca_class))

df = pd.DataFrame(new_dimension_data,columns = ["feature1","feature2","disease"])

df["feature1"] = df["feature1"].astype('float')
df["feature2"] = df["feature2"].astype('float')

#Scatter plot with the new dimensions.
pyplot.title("t-SNE on "+filename)
for i in range(rows):
    pyplot.scatter(df.feature1, df.feature2, c=df.disease)
pyplot.show(block =True)


