eps = 0.0001;
hole_attach_height = 0.02;
outer_radius = 0.25;
hole_radius = 0.045;
thickness = 0.05;
hole_twist = 0;
//other params
inner_radius = outer_radius - thickness;
fn = 12;
fnh = 6;

difference() {
    union() {
        difference() {
            //sphere
            translate([0,0,hole_radius*2 + hole_attach_height])
            sphere(r=outer_radius, $fn = fn);
            
            //shelling
            translate([0,0,hole_radius*2 + hole_attach_height])
            sphere(r=inner_radius, $fn = fn);
            
            translate([0,0,-outer_radius])
            cylinder(h=outer_radius, r= outer_radius, $fn = fn);
        }
        
        //attach segment
        cylinder(h=hole_radius*2 + hole_attach_height, r = outer_radius, $fn = fn);
    }
    
    translate([0,0,-outer_radius -eps/2])
    cylinder(h=hole_radius*2 + hole_attach_height + outer_radius + eps, r = inner_radius, $fn = fn);
    
//    translate([0,0,hole_attach_height + hole_radius]) 
//    rotate([0,0,hole_twist]){
//        rotate([90,0,0])
//        cylinder(h=outer_radius * 2 + eps, r = hole_radius, center=true, $fn = fnh);
//        rotate([0,90,0])
//        cylinder(h=outer_radius * 2 + eps, r = hole_radius, center=true, $fn = fnh);
//    }
    
}