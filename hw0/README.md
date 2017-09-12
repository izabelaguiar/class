
Note that matrix used here to compute first derivatives is a second order accurate
  centered approximation. At the first and last points we implement a second order 
    accurate forward and backward approximation, respectively. 
    All approximations are derived in 
    
        [1] LeVeque, R. "Finite Difference Methods for Ordinary and Partial Differential
            Equations". SIAM (2007), Philadelphia, PA.
            
   The centered approximation is given by eq. (1.3) and the for/back(wards) by eq. (1.11)

Note that matrix used here to compute second derivatives is a second order accurate
    centered approximation. At the first and last points we implement a second order 
    accurate forward and backward approximation, respectively. 
    All approximations are derived in 
    
        [1] LeVeque, R. "Finite Difference Methods for Ordinary and Partial Differential
            Equations". SIAM (2007), Philadelphia, PA.
           
 The centered approximation is given by eq. (1.14) and the for/back(wards) was derived by implementing eq. (1.11) twice.
    
    u'(x) = \frac{(3u(x) - 4u(x+h) + u(x+2h))}{h_i + h_{i+1}}
    u''(x) = \frac{(3u'(x) - 4u'(x+h) + u'(x+2h))}{h_i + h_{i+1}}
            = \frac{1}{h_i+h_{i+1}}(3(\frac{(3u(x) - 4u(x+h) + u(x+2h))}{h_i + h_{i+1}}
                                  -4(\frac{(3u(x+h) - 4u(x+2h) + u(x+3h))}{h_{i+1} + h_{i+2}}
                                  +(\frac{(3u(x+2h) - 4u(x+3h) + u(x+4h))}{h_{i+2} + h_{i+3}}))
           =\frac{1}{h_i+h_{i+1}}(\frac{9}{h_i+h_{i+1}}u_i -(\frac{12}{h_i+h_{i+1}})u_{i+1} 
                                  +(\frac{3}{h_i+h_{i+1}}+\frac{16}{h_{i+1}+h_{i+2}}+\frac{3}{h_{i+2}+h_{i+3}})u_{i+2}
                                  -(\frac{4}{h_{i+1}+h_{i+2}}+\frac{4}{h_{i+2}+h_{i+3}})u_{i+3}
                                  +\frac{1}{h_{i+2}+h_{i+3}}u_{i+4})
                                  
 We then use the test function
      u(x)=sin(x), u'(x)=cos(x), u''(x)=-sin(x)
 to conduct the following numerical experiments.
 
 We first compute the numerical approximation to both the first and second derivative using evenly and unevenly spaced points:
 
 ![alt text](https://github.com/izabelaguiar/class/blob/master/hw0/AS00_num_v_comp.png)
  ![alt text](https://github.com/izabelaguiar/class/blob/master/hw0/AS00_num_v_comp_uneven.png)
 
 We see that the numerical approximations to both of the analytic solutions for the first and second derivative align fairly well
 "in the eyeball norm". This is true for both the evenly spaced points and unevenly spaced points. This observation is not as clear at the endpoints in the second derivative approximation.
 
 To verify that our numerical approximations are indeed second order accurate we conduct a convergence study for decreasing values of h:
 Note that we use the infinity norm to compute the absolute error between the numerical and analytic solutions.
 
  ![alt text](https://github.com/izabelaguiar/class/blob/master/hw0/AS00_cost_v_ac_first.png)
  ![alt text](https://github.com/izabelaguiar/class/blob/master/hw0/AS00_cost_v_ac_second.png)
  
 We see that despite the visibile errors of the endpoints in the first images, we still observe second order convergence in the convergence studies. 
