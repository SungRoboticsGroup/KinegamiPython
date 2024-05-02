//facet number
$fn = 100;

tolerance = 0.01;

//variables
totalLength = 1;
outerRadius = 0.25 - tolerance;
innerRadius = 0.2;
axleRadius = 0.07;

//fixed variables
tolerance = 0.005;
stringHoleRadius = 0.045;
gridHoleSize = 0.03;
thickness = outerRadius - innerRadius;
edgeWidth = outerRadius*2 - 0.15;
topJointBodyLength = 0;
lengthConnect = 0.2;

//grid
module hole(_h) {
    cylinder(h = _h, r = gridHoleSize);
}

// Function to create the circular grid of holes
module circular_grid(_h, _r, _n) {
    for (i = [0:_n-1]) {
        angle = i * 360 / _n;
        x_pos = _r * cos(angle);
        y_pos = _r * sin(angle);
        translate([x_pos, y_pos, -tolerance]) {
            hole(_h);
        }
    }
}

module hole_grid() {
    grid_h = 0.05;
    difference() {
        translate([0,0, stringHoleRadius*4])
        cylinder(h = grid_h, r = outerRadius, center=false);
        translate([0,0,stringHoleRadius*4])
        for (i = [0 : gridHoleSize * 2.5 : outerRadius - gridHoleSize * 2]) {
            circular_grid(grid_h + tolerance*2, i, ceil(6.5*i/(gridHoleSize*3.8)));
        }
    }
}

//topJoint
module jointConnect() {
    difference() {
        union(){
            difference() {
                cylinder(lengthConnect, innerRadius, innerRadius);
                
                translate([0, 0, -0.1])cylinder(1, innerRadius - thickness/2, innerRadius - 
                thickness/2);
            }
        }
        
        //connection screw holes
        union() {
            translate([0, 0.25, lengthConnect - stringHoleRadius*2.5])rotate([90, 90, 0])cylinder(2,
            stringHoleRadius, stringHoleRadius);
        
            translate([-0.25, 0, lengthConnect - stringHoleRadius*2.5])rotate([0, 90, 0])cylinder(2,  
            stringHoleRadius, stringHoleRadius);
        }
        
    }
}

module topJointBody(bodyLength) {
    topJointBodyLength = bodyLength;
    difference() {
        union() {
            difference() {
                union() {
                    difference() {
                        union() {
                            difference() {
                                translate([0, 0, lengthConnect])cylinder(bodyLength, 
                                outerRadius, outerRadius);
                               
                                //string holes
                                translate([0, outerRadius, bodyLength/4 + 
                                lengthConnect])rotate([90, 0, 0])
                                cylinder(outerRadius*2 + 0.25, stringHoleRadius, 
                                stringHoleRadius);
                               
                            }
                        }
                        
                        //cube shelling
                        translate([-edgeWidth/2, -outerRadius, bodyLength/4 + 
                        lengthConnect + 0.1])
                        cube([edgeWidth, 100, 0.5]);
                    }
                }
                //axleHoles
                union() {
                    translate([-outerRadius, 0, bodyLength/2 + (bodyLength/4 + 
                    lengthConnect)]) rotate([0, 90, 0])cylinder(outerRadius + 0.4, 
                    axleRadius, axleRadius);
                    
                    //shelling
                    translate([0, 0, lengthConnect+0.07])cylinder(bodyLength - 0.30, 
                    outerRadius - thickness/2 , outerRadius - thickness/2);
                }
                
            }
        }
        
        //hollow shell
        cylinder(1, innerRadius-thickness/2, innerRadius-thickness/2);     
       
    }
   
    
}

module axle(radius, y) {
    translate([0.5, y, 0])cylinder(outerRadius*2 + 0.1, radius, radius);
}

module bottomJoint(bodyLength) {
    topHalf = edgeWidth;
    bottomHalf = bodyLength - topHalf;
    difference() {
        union() {
            difference() {
                union() {
                    difference() {
                        union() {
                            difference() {
                                union() {
                                    difference() {
                                        translate([1, 0, 0])cylinder(
                                        topJointBodyLength - topJointBodyLength/4 + 
                                        lengthConnect+0.1, edgeWidth/2 - tolerance/2, 
                                        edgeWidth/2 - tolerance/2);
                                        
                                        //axleHoles
                                        translate([1-edgeWidth/2-0.025, 0, 
                                        topJointBodyLength/4 + 0.1])rotate([0, 90, 0 
                                        ])cylinder(outerRadius*2 + 0.25, 
                                        axleRadius, axleRadius);
                                    }
                                    
                                    translate([1, 0, topJointBodyLength - 
                                    topJointBodyLength/4 + lengthConnect+0.1])cylinder
                                    (stringHoleRadius*3 + 0.25, outerRadius, 
                                    outerRadius);
                                }
                                
                                //stringHoles
                                union() {
                                    translate([1, edgeWidth, topJointBodyLength - 
                                    topJointBodyLength/4 + lengthConnect+0.2])rotate([
                                    90, 90, 0])cylinder(edgeWidth+outerRadius, 
                                    stringHoleRadius, stringHoleRadius);
                                    
                                    translate([1, outerRadius, lengthConnect])rotate([90, 90, 0])cylinder(outerRadius*2 + 0.25, stringHoleRadius, stringHoleRadius);
                                    }
                            }
                        }
                    
                        //connecting holes
                        union() {
                            translate([1, outerRadius, topJointBodyLength - 
                            topJointBodyLength/4 + lengthConnect+0.4])rotate([90, 90, 
                            0])cylinder(outerRadius*2,stringHoleRadius, 
                            stringHoleRadius);
                            
                            translate([1 - outerRadius, 0, topJointBodyLength - 
                            topJointBodyLength/4 + lengthConnect+0.4])rotate([0, 90, 
                            0])cylinder(outerRadius*2, stringHoleRadius, 
                            stringHoleRadius);
                        }
                    }
                }
            
            //shellingAll
            translate([1, 0, -0.02])cylinder(bodyLength *2, edgeWidth/2-0.005 - 
            thickness/2, edgeWidth/2-0.005 - thickness/2);
                
            }
        }
        
        //shelling of top (extra)
        translate([1, 0, topJointBodyLength - topJointBodyLength/4 + lengthConnect+0.125])cylinder(bodyLength, innerRadius + 0.025, innerRadius + 0.025);
    }
}



jointConnect();
topJointBody(0.45);
mirror([0, 0, 180])translate([0, 0, -((stringHoleRadius*3 + 0.2) + (topJointBodyLength - topJointBodyLength/4 + lengthConnect+0.15))])bottomJoint(0.45);

translate([1, 0, 0])hole_grid();

//smaller axle
axle(axleRadius-tolerance/0.9, 0.7);
//small axle
axle(axleRadius-tolerance, 0);
//bigger axle
axle(axleRadius - tolerance/2, 0.5);
//exne bigger axle
axle(axleRadius - tolerance/3, -0.5);
