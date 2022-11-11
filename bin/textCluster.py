# from: https://towardsdatascience.com/how-to-easily-cluster-textual-data-in-python-ab27040b07d8
# and: https://towardsdatascience.com/an-approach-for-choosing-number-of-clusters-for-k-means-c28e614ecb2c

import numpy as np
import pickle
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.cluster import KMeans
import sys
    
def cluster_text(text):
    from sklearn.feature_extraction.text import TfidfVectorizer
    vectorizer = TfidfVectorizer(stop_words={'english'})
    X = vectorizer.fit_transform(text)

    Sum_of_squared_distances = []
    K = range(2,6)
    # for k in K:
    #    print(k)
    #    km = KMeans(n_clusters=k, max_iter=200, n_init=10)
    #    km = km.fit(X)
    #    Sum_of_squared_distances.append(km.inertia_)
    
    inertia_o = np.square((X - X.mean(axis=0))).sum()
    for k in K:
       print(k)
       km = KMeans(n_clusters=k, max_iter=200, n_init=10)
       km = km.fit(X)
       
       alpha_k=0.02
    
       Sum_of_squared_distances.append((k,km.inertia_ / inertia_o + alpha_k * k))
    
    results = pd.DataFrame(Sum_of_squared_distances, columns = ['k','Scaled Inertia']).set_index('k')
    best_k = results.idxmin()[0]
    print("best:")
    print(best_k)
#     inertia_o = np.square((scaled_data - scaled_data.mean(axis=0))).sum()
#     # fit k-means
#     kmeans = KMeans(n_clusters=k, random_state=0).fit(scaled_data)
#     scaled_inertia = kmeans.inertia_ / inertia_o + alpha_k * k
    
    plt.plot(results, Sum_of_squared_distances, 'bx-')
    plt.xlabel('k')
    plt.ylabel('Sum_of_squared_distances (adjusted)')
    plt.title('Elbow Method For Optimal k (adjusted)')
    #plt.show()
    plt.savefig("mygraph.png")

    # print('How many clusters do you want to use?')
    # best_k = int(input())
    model = KMeans(n_clusters=best_k, init='k-means++', max_iter=200, n_init=10)
    model.fit(X)
    
    pickle.dump(model, open("model.pkl", "wb"))
    
    # can load later with: 
    model = pickle.load(open("model.pkl", "rb"))

    labels=model.labels_
    clusters=pd.DataFrame(list(zip(text,labels)),columns=['title','cluster'])
    #print(clusters.sort_values(by=['cluster']))

    for i in range(best_k):
        print(clusters[clusters['cluster'] == i])
        
    print("predict:")    
        
    print(model.predict(vectorizer.transform(["My data read from the"])))
    
    return

def main():
    print(sys.argv) # this object will have all of your arguments in it
    print(sys.argv[0]) # this is alwys the name/path? of the python file
    # some_function_in_your_code(sys.argv[0])

    df = pd.read_csv(sys.argv[1], sep=" ")
    print(df)
    
    data = df.iloc[:,2]
    
    print(data)
    
    #data = ["My data read from the Web", "Test here", "test again", "test yet again","Test here", "test again", "test yet again"]

    print(data)

    cluster_text(data)

    #print(modified_data)


if __name__ == "__main__":

    main()
    