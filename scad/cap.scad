eps = 0.0001;
hole_attach_height = 0.02;
inner_radius = 0.2;
outer_radius = 0.25;
hole_radius = 0.045;
//other params
fn = outer_radius * 1000;
fnh = hole_radius * 1000;

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
    
    translate([0,0,hole_attach_height + hole_radius]) {
        rotate([90,0,0])
        cylinder(h=outer_radius * 2 + eps, r = hole_radius, center=true, $fn = fnh);
        rotate([0,90,0])
        cylinder(h=outer_radius * 2 + eps, r = hole_radius, center=true, $fn = fnh);
    }
    
}