**Feedback into the network:**

If we want pattern [1,0,1,0], but we're getting voltage output [0.9,0.7,0.8,0.1] - we can reduce the amount of voltage put into the first electrode, so that the paths for [0.9,0.7] might reduce, leading to reduced voltage of 0.7, maybe to [0.9,0.5].



Then keep doing this until only the target patterns are shown. 

Set a threshold for when the electrode is 'active'. If below threshold, increase voltage by $\delta$V (e.g. 5%). 

Vice versa, if the wrong electrode is 'active', reduce voltage by $\delta$V (e.g. 5%).



<u>**Mac's Ideas:**</u>

- **ADAM - use for gradient descent**:
  - Adaptive Moment Estimate



**Context specificity problem:**

- Have input, output, and a context layer (X or Y) to get two different outputs. E.g. if input + X = one output, but same input + Y = another output

- E.g. X-OR PROBLEM:![image-20210519143030958](C:\Users\61424\AppData\Roaming\Typora\typora-user-images\image-20210519143030958.png)

What are the different logical gates we can implement (AND/OR/NOT/NAND/NOR/XOR)?

- Can we set up a paper that shows we can do this just by changing electrode voltages?



