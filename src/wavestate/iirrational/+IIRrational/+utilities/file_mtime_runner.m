%% @deftypefn  {Function File} {} linux_inotify_wait(dir, func)
%% @deftypefnx {Function File} {} linux_inotify_wait(dir, fname)
%% Wait for any files to change in directory, then run 'func()' or run a script file in the base workspace
%% Requires inotify-tools or similar installed (package search for inotify)
%%
%% @example
%% 
%% @end example
%%
%% @noindent
%% @end deftypefn
function T = file_mtime_runner(directory, func)


  function cleanup()
    disp('DONE DONE DONE');
    delete('timer.update');
  end

  system('touch timer.update');
  T = maketimer(func);

  system('./inotify_touch.sh timer.update &');
  A = onCleanup(@cleanup);
  while true
    start(T);
    wait(T);
  end
end

function T = maketimer(func)
%MUST MUST MUST be separate from the file_mtime_runner since it creates a workspace closure that captures the onCleanup object
%https://stackoverflow.com/questions/22898001/matlab-oncleanup-function-not-executed

  F1 = dir('timer.update');
  function event(obj, event)
    F2 = dir('timer.update');
    if isempty(F2)
      %finish the while loop in the other function ideally
      return
    end
    if F2.datenum > F1.datenum
      F1 = F2;
      if ischar(func)
        evalin('base', ['clear ' func]);
        evalin('base', ['run ' func]);
      else
        func();
      end
    end
  end

  T = timer('Name', 'mtime_watch', ...
            'BusyMode', 'drop', ...
            'ExecutionMode', 'singleShot', ...
            'Period', 2.0, ...
            'StartDelay', 2.0, ...
            'TimerFcn', @event ...
           );
end
