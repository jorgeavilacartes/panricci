# import numpy as np
# import faiss

# import logging

# class Index:

#     def __init__(self, embeddings: np.ndarray):
#         """
#         Args:
#             embeddings (np.ndarray): array with embeddings, one for each row
#         """        
#         # create index
#         self.index = self._create_index(embeddings)

#     def query(self, embeddings, k=10):
#         """return the k nearest neighbors of a set of embeddings

#         Args:
#             embeddings (_type_): a 2-dimensional array with embeddings, one for each row
#             k (int, optional): number of nearest neighbors to return. Defaults to 10.

#         Returns:
#             _type_: tuple of arrays, 
#             the first one contain the distance from closer to distant 
#             the second one contain the id in the index of the closest embeddings 
#             each row correspond to a query embedding, same order as they are input
#         """        
#         distances, indexes = self.index.search(embeddings, k)
#         return distances, indexes
    
#     def _create_index(self,embeddings):
#         "return index built with faiss"
#         dim_embeddings = len(embeddings[0])
#         index = faiss.IndexFlatL2(dim_embeddings)
#         logging.info(f"index created: {index.is_trained}")
#         for emb in embeddings:
#             index.add(emb)
        
#         return index
