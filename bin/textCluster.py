# from: https://towardsdatascience.com/how-to-easily-cluster-textual-data-in-python-ab27040b07d8
# and: https://towardsdatascience.com/an-approach-for-choosing-number-of-clusters-for-k-means-c28e614ecb2c

import numpy as np
import pickle
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.cluster import KMeans
import sys
from sklearn.feature_extraction.text import TfidfVectorizer
from os.path import exists

def cluster_text(text):
    
    vectorizer = TfidfVectorizer(stop_words={'english'})
    X = vectorizer.fit_transform(text)

    Sum_of_squared_distances = []
    K = range(1,6) #text.shape[0])
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
        print(i)
        print(clusters[clusters['cluster'] == i])
        
    print("predict:")    
        
    findID = model.predict(vectorizer.transform(["rccg"]))
    
    print(findID[0])
    
    #print(findID)
    
    cluster = clusters[clusters['cluster'] == findID[0]]
    
    print(type(cluster))
    print(cluster['title'].loc[cluster.index[0]])
    
    
    
    #print(clusters[findID[0]])
    
    
    return

def main():
    
    toSearch=sys.argv[1]
    
    if(exists("model.pkl") and exists("model.data") and exists("vectorizer.pk")):
        # can load later with: 
        model = pickle.load(open("model.pkl", "rb"))
        clusters = pd.read_csv("model.data", sep=",")

        #labels=model.labels_
        #clusters=pd.DataFrame(list(zip(df,labels)),columns=['title','cluster'])
        
        print("predict:")    
        
        vectorizer = pickle.load(open("vectorizer.pk", "rb"))

        findID = model.predict(vectorizer.transform([toSearch]))

        print(findID[0])

        #print(findID)

        cluster = clusters[clusters['cluster'] == findID[0]]

        print(cluster)
        
        print(cluster.iloc[:,2].quantile(q=0.9))
        
        print(cluster.iloc[:,1].loc[cluster.index[0]])
        
        firstCommand = cluster.iloc[:,1].loc[cluster.index[0]]
  
        ratio = len(firstCommand)/len(toSearch)
    
        if(ratio > 0.8 and ratio < 1.25): 
             print(cluster.iloc[:,2].quantile(q=0.9), cluster.iloc[:,3].quantile(q=0.9))
        else: 
            print('noMatch noMatch')
            
    else: 
        vectorizer = TfidfVectorizer(stop_words={'english'})
        print(sys.argv) # this object will have all of your arguments in it
        print(sys.argv[0]) # this is alwys the name/path? of the python file
        # some_function_in_your_code(sys.argv[0])

        df = pd.read_csv("jobRecord.txt", sep=",", header=None)
        #print('df0')
        #print(df)

        df=df[df.iloc[:,2].str.contains('regularSbatch')]

        #print('df')
        #print(df)

        #text = df.iloc[:,0]

        X = vectorizer.fit_transform(df.iloc[:,0])

        Sum_of_squared_distances = []
        K = range(1,6) #text.shape[0])
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
        #with open('vectorizer.pk', 'wb') as fin:
        pickle.dump(vectorizer, open("vectorizer.pk", "wb"))

        #model = pickle.load(open("model.pkl", "rb"))

        labels=model.labels_
        
        print(labels)
        
        print(df.iloc[:,0:4])
        
        clusters=df.iloc[:,[0,12,13]]
        clusters['cluster'] = labels
        
        print(clusters)
        #clusters=pd.DataFrame(list(zip(df.iloc[:,0:4],labels)),columns=['title','mem','time','cluster'])
        
        clusters.to_csv("model.data")
        
        
        for i in range(best_k):
            print(i)
            print(clusters[clusters['cluster'] == i])

        print("predict:")    

        findID = model.predict(vectorizer.transform([toSearch]))

        print(findID[0])

        #print(findID)

        cluster = clusters[clusters['cluster'] == findID[0]]
        print(cluster)
        
        print(cluster.iloc[:,1].quantile(q=0.9))
        
        
            
  
        
        print(type(cluster))
        print(cluster.iloc[:,0].loc[cluster.index[0]])

        firstCommand = cluster.iloc[:,0].loc[cluster.index[0]]
  
        ratio = len(firstCommand)/len(toSearch)
        if(ratio > 0.8 and ratio < 1.25): 
             print(cluster.iloc[:,1].quantile(q=0.9), cluster.iloc[:,2].quantile(q=0.9))
        else: 
            print('noMatch noMatch')
            
    return
    
if __name__ == "__main__":

    main()
    