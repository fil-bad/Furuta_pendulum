function [values,isterminal,direction] = time_event(t,x,tstart)
 %  Don't let t cross zero...a dummy "event" to illustrate how 
 %  one might handle other events in conjunction with the time
 %  constraint.  Not necessary, but I put it in just in case.
 values(1) = t;
 
 %  Don't let integration go for more than 2 seconds.
 values(2) = toc(tstart) < 2;
 isterminal = true(size(values));
 direction = zeros(size(values));
end