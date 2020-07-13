function y = str2syms(symstr)
%STR2SYMS evaluates symstr where symstr is a string 
% representing a symbolic expression.
% Additionally to str2sym, str2syms declares unknown variables
% as symbolic in the caller workspace.
%
% Example:
% syms x y
% y = x;
% z = str2syms('a + y')
% creating symbolic variable a
% z =
% a + x

% very simple tokenizer to find variables
vars = strread(symstr,'%s','delimiter',',:;= ()^*/+-.0123456789@');
% collect all known variables
W = evalin('caller','whos'); %or 'base'
for i = 1:numel(vars)
    v = vars{i};
    if ~isempty(v)
        % create symbolic variable if not already known
        if ~ismember(v,[W(:).name]) & ~exist(v)
        %if ~exist(v)
            disp(['creating symbolic variable ', v])
            assignin('caller', v, sym(v));
        end
    end
end
% convert (single) = to ==
ieq = strfind(symstr,'=');
if ieq
    symstr = [symstr(1:ieq),'=',symstr(ieq+1:end)];
end
% evaluate
y = evalin('caller', symstr);
