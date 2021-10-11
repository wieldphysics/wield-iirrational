

ds=DepStruct();
ds.add_generator('test', @(v)1:10);
ds.add_generator('test2', @(v)2*v.test);
disp(ds.test2)
ds.clear('test')
ds.add_generator('test', @(v)10:20);
disp(ds.test2)
%assert(1);

ds=DepStruct();
ds.test = 1:5;
