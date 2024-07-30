# Three_Body_Problem_Code
This repository contains the code used for my summer project on the three-body problem, which resulted in the plots published on [my reserach website](https://henry-yip.github.io/tagsbycontent/#research).
References and descriptions are embedded within the code itself. For any inquiries, please do not hesitate to email me.

# Important Convention
The xyz files required in most of the provided codes require the below format: 

```plaintext
3
Point = 1
Central_Planet    -6.942220069422201e-05 -0.00010654559946727201 0.0
First_Planet      1.0000347048503473 5.3272799733677604e-05 0.0
Second_Planet     -0.9999652826496529 5.32727997335944e-05 0.0
3
Point = 2
Central_Planet    -0.00013884440555377584 -0.00021309119573817583 0.0
First_Planet      1.000069397202777 0.00010654559786975359 0.0
Second_Planet     -0.999930552797223 0.00010654559786842221 0.0
3
Point = 3
Central_Planet    -0.0002082666187439927 -0.0003196367856163424 0.0
First_Planet      1.0001040770593717 0.0001598183928115413 0.0
Second_Planet     -0.9998958104406274 0.00015981839280480113 0.0
```

 - The initial 3 indicates the number of particles (always 3 in a three-body problem).
 - Each "Point" refers to a specific timestep.
 - The subsequent lines provide the coordinates corresponding to the positions of the particles at each timestep.
