{
    // element by element vector arithmatic
    vadd(v1,v2) :: std.map(function(i) v1[i]+v2[i],std.range(0,std.length(v1)-1)),
    vmul(v1,v2) :: std.map(function(i) v1[i]*v2[i],std.range(0,std.length(v1)-1)),
    vsub(v1,v2) :: std.map(function(i) v1[i]-v2[i],std.range(0,std.length(v1)-1)),

    // return the direction of v
    vdir(v) :: $.scale(v, 1.0 / $.mag(v)),

    // Return a vector displaced along direction of d by length l from v.
    vshift(v,d,l) :: $.vadd(v, $.scale($.vdir(d), l)),

    // Return new vector of absolute value of input elements 
    vabs(v) :: std.map(std.abs, v),

    // Reduce vector by summing elements.
    sum(v) :: std.foldl(function(n,x) n+x, v, 0),

    // Magnitude/length of a vector 
    mag(v) :: std.sqrt(self.sum(std.map(function(x) x*x, v))),

    // Reduce vector by multiplication.
    mul(v) :: std.foldl(function(n,x) n*x, v, 1),

    // Multiply a scalar to a vector, element by element
    scale(v,s) :: std.map(function(e) e*s, v),

    // Add a scalar to a vector, element by element
    increase(v,s) :: std.map(function(e) e+s, v),

    // Return 3-vector in "3D point" format
    topoint(v) :: { x:v[0], y:v[1], z:v[2] },

    // return 3-vector given a "3D point" 
    frompoint(p) :: [ p.x, p.y, p.z ],

    // Calculate the volume of a rectangular solid defined by two
    // points giving extreme corners.
    volume(v1, v2) :: $.mul($.vabs($.vsub(v1,v2))),
}
