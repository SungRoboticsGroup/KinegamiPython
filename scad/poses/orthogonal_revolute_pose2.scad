tolerance = 0.015;
hole_radius = 0.045;
grid_hole_radius = 0.03;
outer_radius = 0.3;
thickness = 0.05;
hole_attach_height = 0.02;
attach_thickness = 0.02;
bottom_rotation_height = 1;
top_rotation_height = 0.7;
next_hole_radius = 0.045;
next_hole_attach_height = 0.02;
next_inner = 0.2;
max_turn = 60;
min_turn = -30;
hole_twist = 0;

//parameters that don't have to be defined explicitly
eps = tolerance/100;
inner_radius = outer_radius - thickness;
axle_offset = attach_thickness;
bottom_joint_attach_radius = outer_radius - hole_attach_height - grid_hole_radius*2 - thickness;
fn = 12;
fnh = 6;

//solve for b = axle height
max_turn_b = outer_radius / tan((180 - max_turn)/2);
min_turn_b = outer_radius / tan((180 - min_turn)/2);
axle_radius = bottom_joint_attach_radius/3;

grid_hole_height = hole_radius*2 + hole_attach_height * 2;
bottom_joint_grid_clearance = bottom_rotation_height - (hole_attach_height*2 + hole_radius*2 + grid_hole_height + max_turn_b + attach_thickness);

module angled_cylinder(angle_deg_y, _h, _r, angle_deg_x = 0) {
    rotate([angle_deg_x, angle_deg_y, 0]) { // Rotate the cylinder 
            cylinder(h = _h, r = _r, center = true, $fn = fnh);
        }
}

module joint_top() {
    operating_height = hole_attach_height*2 + attach_thickness + hole_radius*2 + grid_hole_height + bottom_joint_grid_clearance + max_turn_b;
    
    top_attach_height = axle_offset + axle_radius + max_turn_b;
    cyl_height = top_rotation_height - max_turn_b;

    difference() {
        translate([0,0, operating_height]) 
        difference() {
            union() {
                translate([0,0,- axle_radius - axle_offset + top_attach_height/2]) 
                {
                    
                    //attach points
                    axle_attach_thickness = outer_radius - thickness/2 - bottom_joint_attach_radius - tolerance ;
                    
                    attach_distance = bottom_joint_attach_radius + axle_attach_thickness/2 + tolerance;

                    translate([attach_distance, 0, 0]) {
                        translate([0,0,-top_attach_height/2 + axle_radius + axle_offset])
                        rotate([0,90,0])
                        cylinder(h=axle_attach_thickness,r=axle_radius + axle_offset,center=true, $fn = fn);
                        translate([0,0,(axle_radius + axle_offset)/2])
                    cube([axle_attach_thickness, (axle_radius + axle_offset)*2, top_attach_height - axle_offset - axle_radius + attach_thickness],center=true);
                    }
                    
                    translate([-attach_distance, 0, 0]) {
                        translate([0,0,-top_attach_height/2 + axle_radius + axle_offset])
                        rotate([0,90,0])
                        cylinder(h=axle_attach_thickness,r=axle_radius + axle_offset,center=true, $fn = fn);
                        translate([0,0,(axle_radius + axle_offset)/2])
                    cube([axle_attach_thickness, (axle_radius + axle_offset)*2, top_attach_height - axle_offset - axle_radius + attach_thickness],center=true);
                    }  
                    
                    //outer cylinder
                    translate([0,0,top_attach_height/2]) {
                        difference() {
                            cylinder(h=cyl_height,r=outer_radius,$fn=fn);
                            translate([0,0,-eps/2])
                            cylinder(h=cyl_height + eps,r=inner_radius,$fn=fn);
                        }
                        
                        //enforce min turn
                        translate([0,0,-(max_turn_b - min_turn_b)])
                        difference() {
                            cylinder(h=max_turn_b - min_turn_b,r=outer_radius,$fn=fn);
                            
                            translate([-outer_radius,-outer_radius,-eps/2])
                            cube([outer_radius*2,outer_radius - axle_radius - axle_offset + eps/2,max_turn_b - min_turn_b + eps]);
                            
                            translate([0,0,-eps/2])
                            cylinder(h=max_turn_b - min_turn_b + eps,r=inner_radius,$fn=fn);
                       }
                    }       
                } 
                
                
                //attachment segment
                translate([0,0,cyl_height + max_turn_b])
                difference() {
                        union() {
                            //attachment for smaller cylinder
                            translate([0,0,-attach_thickness])
                            cylinder(h=attach_thickness, r=outer_radius, $fn=fn);
                            translate([0,0,-attach_thickness*2])
                            cylinder(h=attach_thickness*2, r1=outer_radius, r2=next_inner - tolerance - thickness, $fn=fn);
                            
                            //small cylinder
                            cylinder(h=next_hole_radius*2 + next_hole_attach_height * 2, r = next_inner - tolerance, $fn=fn);
                        }
                        
                        //shelling
                        translate([0,0, -attach_thickness*2 - eps/2])
                        cylinder(h=next_hole_radius*2 + next_hole_attach_height * 2 + attach_thickness*2 + eps, r = min(inner_radius, next_inner - tolerance - thickness), $fn=fn);
                        
                        //holes for screws
//                        translate([0,0,next_hole_attach_height + next_hole_radius]) {
//                            angled_cylinder(0, (next_inner + thickness + eps) * 2,  next_hole_radius, 90);
//                            angled_cylinder(90, (next_inner + thickness + eps) * 2, next_hole_radius, 0);
//                        }
                    }

            }
            
            //holes for string routing
//             translate([0,-inner_radius + tolerance,max_turn_b]) 
//            rotate([90,0,0])
//            {
//                cylinder(h=thickness + inner_radius + tolerance,r=grid_hole_radius,$fn=fn);
//                translate([0,hole_attach_height + grid_hole_radius*2,0])
//                cylinder(h=thickness + inner_radius + tolerance,r=grid_hole_radius,$fn=fn);
//
//            }
//
//            translate([0,outer_radius + tolerance,min_turn_b])
//            rotate([90,0,0])
//            {
//                cylinder(h=thickness + inner_radius + tolerance,r=grid_hole_radius,$fn=fn);
//                translate([0,hole_attach_height + grid_hole_radius*2,0])
//                cylinder(h=thickness + inner_radius + tolerance,r=grid_hole_radius,$fn=fn);
//
//            }
             

        }

        //axle hole
        translate([0,0,operating_height])
        rotate([0,90,0])
        cylinder(h=outer_radius*2,r=axle_radius,center=true,$fn=fn);
        
    }
}

translate([0,0,-(bottom_rotation_height + top_rotation_height)])
joint_top();