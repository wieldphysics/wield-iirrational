%% @deftypefn  {Function File} {} polynomial ()
%% @deftypefnx {Function File} {} polynomial (@var{a})
%% Create a polynomial object representing the polynomial
%%
%% @example
%% @end example
%%
%% @noindent
%% @end deftypefn
classdef DepStub < handle
  properties
    other
  end
  methods (Access = public)
    function self = DepStub(other)
      self.other = other;
    end

    function [varargout] = subsref(self, ref)
      [varargout{1:nargout}] = subsref(self.other, ref);
    end

    function self = subsasgn(self, ref, val)
      subsasgn(self.other, ref, val);
    end
  end
end
%!test
%! ds=DepStruct();
%! ds.generator('test', @(v)1:10);
%! ds.('test');
%! assert(1);
%! ds.do_test()

