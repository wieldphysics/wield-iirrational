%% @deftypefn  {Function File} {} polynomial ()
%% @deftypefnx {Function File} {} polynomial (@var{a})
%% Create a polynomial object representing the polynomial
%%
%% @example
%% @end example
%%
%% @noindent
%% @end deftypefn
function names = fieldnames(self)
  names = [
          fieldnames(self.values),
          fieldnames(self.generators),
          fieldnames(self.setters),
          methods(self),
          builtin('properties', self),
  ];
  names = unique(names);
end
