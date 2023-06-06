# ECE 208 Homework 9: Bootstrapping and model selection
In this homework assignment, we will use bootstrapping to choose between alternative models and to select among a set of hypotheses trees. This assignment will heavily use your implementation of the [Felsenstein's tree-pruning algorithm](https://en.wikipedia.org/wiki/Felsenstein%27s_tree-pruning_algorithm)<sup>1</sup> for computing the (log-)likelihood under the [General Time Reversible (GTR) model](https://en.wikipedia.org/wiki/Models_of_DNA_evolution#GTR_model_(Tavar%C3%A9_1986))<sup>2</sup> given a tree *T* with known branch lengths and topology. 


## Problem
Given are:

- an alignment, 
- a *main* tree ***T***
- a *main* rate matrix ***R***
- One of the following two:
	- an *alternative* tree ***T<sub>A</sub>***
	- an *alternative* rate matrix ***R<sub>A</sub>***

Assume the main tree/matrix pair has a better log-likelihood than alternative trees and alternative matrices. We want to test whether the difference in the likelihood of the main tree and the alternative tree is significant or not (ditto for the rate matrix). Thus, we want to decide whether we can reject the alternative trees and matrices; if so, we can state that they are inferior to the main tree and/or matrix in a statistically rigorous way. 

You will approach this problem using bootstrapping. 

1. Using bootstrapping, you should compute a distribution over the likelihood scores for the *main* ***(T,R)*** pair. The resulting distribution is your "null" distribution.  
	-  To compute likelihood you can use your code from HW7, but see hints below.
2. Then, you compute the likelihood of the alternative tree ***T<sub>A</sub>*** or the alternative rate matrix ***R<sub>A</sub>***, whichever is given. Call this likelihood ***L<sub>A</sub>***.
3. You will report the probability of observing a likelihood as low or lower than ***L<sub>A</sub>***. You compute this probability with reference to the null distribution computed in step 1. Let's call this the p-value for this input. 

**Hints**:

1. Instead of bootstrapping the alignment, see if you can find a more clever approach that saves time. We leave it to you to figure out what that means. 
2. You are allowed to change your likelihood code from HW7. 
3. If you do this correctly, bootstrapping can be done in one line of code (after a small change to the code from HW7). The entire program could be less than 10 lines of code.
4. Feel free to use numpy functions for random selection. 	    

 
### Input

In this assignment, there is no starter code. You are on your own. 

Provided inputs are:

- A directory called **trees** with a set of 27 alternative trees
- A directory called **params** with a set of 5 alternative rate matrices
- A directory called **example** with 
	- 10 input alignments called ***01.fas -- 10.fas***
	- The main tree called **tree.nwk**
	- The main rate matrix called **gtr_params.txt**

Expected output is:

- A .csv file with exactly 330 lines. Each line will look like this:
	```
	2 TN93_params.txt tree.nwk 0.000100
	``` 
	- The first number is the alignment number, the second field is the name of the parameter file, the third field is the name of the tree file, and the last field is the p-value.
- Each 33 lines correspond to one of the 10 input alignment files. For each alignment file, you report:
	- 5 lines: the results of testing each of the parameter files in the params directory (but using the main tree)
	- 27 lines: the results of testing each of the tree files in the trees directory (but using the main parameter file)
	- 1 line: the result of using the main tree and the main parameter file (as if they were the alternative alignment and tree)



We have provided an [example outputfile](output-example.csv). All you really have to do is to replace the last column of this file with the p-value you computed (without changing the other columns). It's fine if the order of lines in the file is changed. 

## Homework Deliverables
* Your code for implementing the test. We won't run the code. Just provide it for reference. Give us the name of the file and its usage in the writeup. 

* The **output-example.csv** file with values on the last column replaced with results of your calculations. 

* **[writeup.pdf](writeup.pdf):** A write-up that includes the following pieces of information:
    1. Summarize (perhaps using visualization) which of the 27 trees in your judgment are likely to have generated the data. Which can be safely rejected? Which are ambiguous? 
    2. Summarize (perhaps using visualization) which of the 5 models are overall as good as GTR
    3. Which model would *you* choose? Why?
    4. How did you perform bootstrapping? 
    5. What you learned from the exercise (a short paragraph is enough). 
        - (Optional) Feel free to comment briefly on the procedure we used and if you think it has any shortcomings. Are there ways to improve it?

    The whole writeup can be as short as half a page. If you need more space, use it. But don't feel it needs to be long. 

## Grade Breakdown (100 Points)
* **[writeup.pdf](writeup.pdf): 45 Points**
    * *A well-reasoned discussion of what models could have generated the data (i,ii):* 20 points
    * *A well-reasoned selection of the model (iii):* 10 points
    * *An efficient bootstrapping procedure (iv)*: 10 points
    * *Lessons learned (v)*: 5 points. 
        - 10 extra points for an intelligent critique of the approach. 

* **[output-example.csv](output-example.csv): 55 Points**
    * Since bootstrapping is random, we won't expect to have identical results between your runs and ours. But, statistically, your p-values should be close to ours. We will judge the accuracy of your calculations by comparing your results to ours, allowing for random noise. Each of the 330 points will have a 1/6 of a point.

## References
1. Felsenstein J. Maximum Likelihood and Minimum-Steps Methods for Estimating Evolutionary Trees from Data on Discrete Characters. *Systematic Biology*, 22(3):240–249, 1973.
2. Tavaré S. Some probabilistic and statistical problems in the analysis of DNA sequences. *Lectures on Mathematics in the Life Sciences*, 17:57–86, 1986.
