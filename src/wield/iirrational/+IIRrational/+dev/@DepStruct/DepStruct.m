%% @deftypefn  {Function File} {} polynomial ()
%% @deftypefnx {Function File} {} polynomial (@var{a})
%% Create a polynomial object representing the polynomial
%%
%% @example
%% @end example
%%
%% @noindent
%% @end deftypefn
classdef DepStruct < handle
  properties
    values
    values_computed
    values_mark
    generators
    setters
    value_dependencies
    current_func
    current_autodep
    current_mark
    autodeps_setters
    autodeps_generators
    my
  end
  methods (Access = public)
    function self = DepStruct(varargin)
      if nargin > 0 && isa(varargin{1}, class(self))
        copy = p.Results.copy;
      else
        copy = [];
      end

      if not(isempty(copy))
        %has copy constructor
        p = inputParser();
        %need to allow other input
        p.addOptional('copy', []);
        p.addParamValue('lightweight', 1);
        p.gaineepUnmatched = true;
        p.parse(varargin{:});
      end

      self.current_mark = 0;
      if isempty(copy)
        self.values              = struct();
        self.values_computed     = struct();
        self.values_mark         = struct();
        self.value_dependencies  = struct();

        self.generators          = struct();
        self.setters             = struct();

        self.autodeps_setters    = struct();
        self.autodeps_generators = struct();

      else
        self.current_mark = copy.current_mark;

        if p.Results.lightweight
          self.values              = struct();
          self.values_computed     = struct();
          self.values_mark         = struct();
          self.value_dependencies  = struct();
          for f = fieldnames(copy.values)
            f = f{1};
            if copy.values_computed.(f);
              self.values.(f) = copy.values.(f);
              self.values_computed.(f) = copy.values_computed.(f);
              self.values_mark.(f)     = copy.values_mark.(f);
              %TODO, this isn't quite right, some deps should be purged
              self.value_dependencies.(f) = copy.value_dependencies.(f);
            end
          end
        else
          self.values              = copy.values;
          self.values_computed     = copy.values_computed;
          self.values_mark         = copy.values_mark;
          self.value_dependencies  = copy.value_dependencies;
        end

        self.generators          = copy.generators;
        self.setters             = copy.setters;

        self.autodeps_setters    = copy.autodeps_setters;
        self.autodeps_generators = copy.autodeps_generators;

      end

      self.current_func        = [];
      self.current_autodep     = 0;

      self.my                  = DepStub(self);

      if isempty(copy)
        self.build(varargin{:})
      end
    end

    function build(self, varargin)
    end

    function val = compute(self, name)
      %make this protected
      prev_cf = self.current_func;
      prev_adcf = self.current_autodep;
      self.current_autodep = self.autodeps_generators.(name);
      self.current_func = name;
      gen_func = self.generators.(name);

      %send in self.my rather than self since in subclasses subsref will not be working even in anonymous
      %local functions
      newval = gen_func(self.my);
      self.values.(name)          = newval;
      self.values_computed.(name) = 1;
      self.values_mark.(name)     = self.current_mark;
      self.current_func           = prev_cf;
      self.current_autodep        = prev_adcf;
      val = newval;
      %disp(self.value_dependencies)
    end

    function newval = assign(self, name, val)
      prev_cf = self.current_func;
      prev_adcf = self.current_autodep;
      self.current_autodep = self.autodeps_setters.(name);
      self.current_func = name;
      gen_func = self.setters.(name);

      %send in self.my rather than self since in subclasses subsref will not be working even in anonymous
      %local functions
      newval = gen_func(self.my, val);
      self.values.(name)          = newval;
      self.values_computed.(name) = 0;
      self.values_mark.(name)     = self.current_mark;
      self.current_func           = prev_cf;
      self.current_autodep        = prev_adcf;
    end

    function set_raw(self, name, val)
      mark = self.current_mark;
      self.current_mark = self.current_mark + 1;
      newval = assign(self, name, val);
      self.clear_dependent(name, mark);
    end

    function val = get_raw(self, name)
      if isfield(self.values, name)
        val = self.values.(name);
        return
      end
      val = compute(self, name);
    end

    function [varargout] = subsref(self, ref)
      switch ref(1).type
          case '.'
            name = ref(1).subs;

            if self.current_autodep && not(isempty(self.current_func))
              dependencies_for(self, self.current_func, name);
            end

            if isfield(self.values, name)
              newval = self.values.(name);
              if not(isempty(ref(2:end)))
                [varargout{1:nargout}] = subsref(newval, ref(2:end));
              else
                varargout{1} = newval;
              end
              return
            elseif not(isfield(self.generators, name))
              [varargout{1:nargout}] = builtin('subsref', self, ref);
              return
            end

            newval = compute(self, name);
            if not(isempty(ref(2:end)))
              [varargout{1:nargout}] = subsref(newval, ref(2:end));
            else
              varargout{1} = newval;
            end

            return
          otherwise
            error('only indexes as a struct/object');
       end
    end

    function self = subsasgn(self, ref, val)
      switch ref(1).type
        case '.'
          name = ref(1).subs;
          if not(isfield(self.setters, name))
            self.add_setter(name);
          end

          mark = self.current_mark;
          self.current_mark = self.current_mark + 1;
          newval = assign(self, name, val);
          self.clear_dependent(name, mark);

          %this is likely not the logic that we want
          if not(isempty(ref(2:end)))
            subsasgn(newval, ref(2:end), newval);
          end
          return
        otherwise
          error('only indexes as a struct/object');
      end
    end

    function clear(self, name)
      if isfield(self.values, name)
        self.values = rmfield(self.values, name);
        self.values_computed = rmfield(self.values_computed, name);
      end
      self.clear_marked(name, self.current_mark);
    end

    function clear_marked(self, name, mark)
      if isfield(self.values, name)
        self.values = rmfield(self.values, name);
        self.values_computed = rmfield(self.values_computed, name);
      end

      if not(isfield(self.value_dependencies, name))
        return
      end

      deps = self.value_dependencies.(name);
      self.value_dependencies = rmfield(self.value_dependencies, name);
      for dep = fieldnames(deps)
        if self.values_mark.(dep{1}) <= mark
          clear_marked(self, dep{1}, mark);
        end
      end
    end

    function clear_dependent(self, name, mark)
      %make protected
      if not(isfield(self.value_dependencies, name))
        return
      end
      deps = self.value_dependencies.(name);
      self.value_dependencies = rmfield(self.value_dependencies, name);
      for dep = fieldnames(deps)
        if self.values_mark.(dep{1}) <= mark
          clear_marked(self, dep{1}, mark);
        end
      end
    end

    function add_setter(self, varargin)
      p = inputParser();
      %must add validator for it to take a string
      p.addOptional('name', [], @(v) 1);
      p.addOptional('func', []);
      p.addOptional('autodeps', 1);
      p.addOptional('clear', 1);
      p.parse(varargin{:});
      name = p.Results.name;
      func = p.Results.func;

      if isempty(name)
        name = strsplit(func2str(func), '/');
        name = name{end};
      elseif isempty(func)
        %ok, then we were only given a name
        if isa(name, 'function_handle')
          func = name;
          name = strsplit(func2str(func), '/');
          name = name{end};
        else
          %if name is a string, then we just assign a simple setter
          func = @(ds, val) val;
        end
      end
      self.autodeps_setters.(name) = p.Results.autodeps;
      self.setters.(name) = func;
      if p.Results.clear
        self.clear(name);
      end
    end

    function add_generator(self, varargin)
      p = inputParser();
      %must add validator for it to take a string
      p.addOptional('name', [], @(v) 1);
      p.addOptional('func', []);
      p.addOptional('autodeps', 1);
      p.addOptional('clear', 1);
      p.parse(varargin{:});
      name = p.Results.name;
      func = p.Results.func;

      if isempty(func)
        assert(isa(name, 'function_handle'));
        func = name;
        name = strsplit(func2str(func), '/');
        name = name{end};
      elseif isempty(name)
        name = strsplit(func2str(func), '/');
        name = name{end};
      end
      self.autodeps_generators.(name) = p.Results.autodeps;
      self.generators.(name) = func;
      if p.Results.clear
        self.clear(name);
      end
    end

    function dependencies(self, varargin)
      if isempty(self.current_func)
        error('dependencies only callable from inside a setter or generator (use dependencies_for perhaps)');
      end
      dependencies_for(self, self.current_func, varargin{:});
    end

    function dependencies_for(self, name, varargin)
      for value = varargin
        self.value_dependencies.(value{1}).(name) = 1;
      end
      %disp(self.value_dependencies)
    end
  end
end

%!test
%! ds=DepStruct();
%! ds.add_generator('test', @(v)1:10);
%! ds.test;
%! assert(ds.test == 1:10);

%!test
%! ds=DepStruct();
%! ds.add_generator('test', @(v)1:10);
%! ds.add_generator('test2', @(v)2*v.test);
%! assert(ds.test2 == (2*(1:10)));
%%! disp(ds.test2);
%! ds.clear('test')
%! ds.add_generator('test', @(v)10:20);
%! assert(ds.test2 == (2*(10:20)));
%! assert(1);

%!test
%! ds=DepStruct();
%! ds.add_setter('test', @(v, val) val);
%! ds.test = 1:5;
%! assert(ds.test == (1*(1:5)));

