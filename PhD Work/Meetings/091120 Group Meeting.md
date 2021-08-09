- Rerun case studies with binary model and compare accuracy

![image-20201109171041679](C:\Users\61424\AppData\Roaming\Typora\typora-user-images\image-20201109171041679.png)

- RURU conductance vs voltage idea:
  - Rerun with T = 10 instead of T = 5:

```python
kev = maxConductance[case][:4000].reshape(2,2000)
np.mean(np.std(kev, axis = 0))
```

-ask joel ideas re:tunnelling 