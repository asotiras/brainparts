## Brainlets

NMF Code for calculating components from Brain imaging data.

The first matlab functions (opnmfm and opnmf_mem.m) implement the orthogonal projective non-negative matrix factorization. The only difference between them is that the order of the matrix multiplications in the multiplicative update rule of opnmf_mem.m function has been rearranged so that high dimensional imaging data can used. For low dimensional data, the function opnmf.m should be faster. Note that in both functions, at every iteration, really low values are thresholded out as they result in increased computational cost.

NNSVD.m initializes the factorization. In our experiments, we have found it to be much more efficient than random initialization and the first (default) option is deterministic. It works by operating a dual SVD decomposition. The SVD may require a lot of memory resources depending on the size and number of images. Thus, a fast randomized SVD procedure (flag=4) is often employed. The randomized procedure makes use of randomized principal component analysis (randpca.m file).  In our experiments, we have found the approximation to be quite accurate and useful for practical purposes since it decreases the overall memory requirements of the code.

The default values are the ones we ypically use in our experiments. Note that we let the algorithm run for many iterations to ensure convergence. As a consequence, the algorithm may run for a long time (especially for high-dimensional imaging data). To counter possible failures, wesave intermediate results that can be used to restart the experiment in case something goes wrong.

----

## Reference

Sotiras, Aristeidis, Susan M. Resnick, and Christos Davatzikos. **Finding imaging patterns of structural covariance via non-negative matrix factorization.** *NeuroImage* 108 (2015): 1-16.

Yang, Zhirong, and Erkki Oja. **Linear and nonlinear projective nonnegative matrix factorization.** *IEEE Transactions on Neural Networks* 21, no. 5 (2010): 734-749.

Boutsidis, Christos, and Efstratios Gallopoulos. **SVD based initialization: A head start for nonnegative matrix factorization.** *Pattern Recognition* 41, no. 4 (2008): 1350-1362.

Halko, Nathan, Per-Gunnar Martinsson, and Joel A. Tropp. **Finding structure with randomness: Probabilistic algorithms for constructing approximate matrix decompositions.** *SIAM review* 53, no. 2 (2011): 217-288.

----

## Papers using this code to identify patterns in neuroimaging data

Sotiras, Aristeidis, Jon B. Toledo, Raquel E. Gur, Ruben C. Gur, Theodore D. Satterthwaite, and Christos Davatzikos. **Patterns of coordinated cortical remodeling during adolescence and their associations with functional specialization and evolutionary expansion.** *Proceedings of the National Academy of Sciences* 114, no. 13 (2017): 3527-3532.

Marieta Pehlivanova, Daniel H. Wolf, Aristeidis Sotiras, Antonia Kaczkurkin, Tyler M. Moore, Rastko Ciric, Philip A. Cook, Angel Garcia de La Garza, Adon Rosen, Kosha Ruparel, Anup Sharma, Russell T. Shinohara, David R. Roalf, Ruben C. Gur, Christos Davatzikos, Raquel E. Gur, Joseph W. Kable and Theodore D. Satterthwaite. **Diminished cortical thickness is associated with impulsive choice in adolescence.** *Journal of neuroscience* (2018): 2200-17.

Nassar, Rula, Antonia N. Kaczkurkin, Cedric Huchuan Xia, Aristeidis Sotiras, Marieta Pehlivanova, Tyler M. Moore, Angel Garcia de La Garza, David R. Roalf, Adon F. G. Rosen, Scott A. Lorch, Kosha Ruparel, Russell T. Shinohara, Christos Davatzikos, Ruben C. Gur, Raquel E. Gur, and Theodore D. Satterthwaite. **Gestational Age is Dimensionally Associated with Structural Brain Network Abnormalities Across Development.** *Cerebral Cortex* (2018).

Varikuti, Deepthi P., Sarah Genon, Aristeidis Sotiras, Holger Schwender, Felix Hoffstaedter, Kaustubh R. Patil, Christiane Jockwitz, Svenja Caspers, Susanne Moebus, Katrin Amunts, Christos Davatzikos, and Simon Eickhoff. **Evaluation of non-negative matrix factorization of grey matter in age prediction.** *NeuroImage* 173 (2018): 394-410.

Antonia Kaczkurkin, Rula Nassar, Cedric Xia, Aristeidis Sotiras, Marieta Pehlivanova, Tyler M. Moore, Angel Garcia de la Garza, David R. Roalf, Adon Rosen, Scott Lorch, Kosha Ruparel, Russell T. Shinohara, Christos Davatzikos, Ruben C. Gur, Raquel E. Gur, Theodore D. Satterthwaite. **A Dimensional Measure of Prematurity is Associated With Structural Brain Network Abnormalities in Children, Adolescents, and Young Adults.** *Neuropsychompharmachology* 43 (2017): S304-S305

Yin Chen, Aristeidis Sotiras, Ilya Nasrallah, Rizwan Akhtar, Jacqueline Rick, Alice Chen-Plotkin, John Trojanowski, Daniel Weintraub, Christos Davatzikos, Jacob Dubroff. **Non-negative matrix factorization evaluation of patterns of brain Aβ deposition in Parkinson’s disease, Alzheimer’s disease and normal controls on \[<sup>18F<sup>\]florbetapir PET.** *Journal of Nuclear Medicine* 59 (2018):484-484
  
----
