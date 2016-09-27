%% gaussify(A,b)
% Gaussian elimination on matrix A
% A = matrix to be eliminated,
% b = right hand side (defaults to eye(size(A)).

function gaussify(A,b)
A = sym(A);
if nargin < 2 || numel(b) < 2
    B = eye(size(A));
else
    B = b;
end
B = sym(B);

str = '';

[r,~] = size(A);

ch = true;
% Downward Elimination
for j = 1:r
    if ch
        gms = gmprint(A,B);
        str = [str sprintf('&\\begin{gmatrix}[p]\n') gms];
    end
    if j < r
        if ch
            str = [str sprintf('\\rowops\n')];
        end
        ch = false;
        for k=j+1:r
            multip = -A(k,j)/A(j,j);
            
            if multip ~= 0
               ch = true;
            else
                continue
            end
            
            A(k,:) = multip*A(j,:) + A(k,:);
            B(k,:) = multip*B(j,:) + B(k,:);
            
            str = [str sprintf('\\mult{%d}{%s}\n\\add{%d}{%d}\n', j-1,['\cdot ' sprintf(format_mult(multip))],j-1,k-1)];
        end
        if ch
            str = [str sprintf('\\end{gmatrix}\n\\\\ \\implies \n')];
        end
    else
        if ch
            str = [str sprintf('\\rowops\n')];
        end
        ch = false;
        for k = r-1:-1:1
            multip = -A(k,j)/A(j,j);
            
            if multip ~= 0
                ch = true;
            else
                continue
            end
            
            A(k,:) = multip*A(j,:) + A(k,:);
            B(k,:) = multip*B(j,:) + B(k,:);
            
            str = [str ...
                sprintf('\\mult{%d}{%s}\n\\add{%d}{%d}\n', ...
                j-1,['\cdot ' sprintf(format_mult(multip))],j-1,k-1)];
        end
        str = [str sprintf('\\end{gmatrix}\n')];
        if ch
            str = [str sprintf('\\\\\\implies')];
        end
    end
end

% Upwards
for j = r-1:-1:1
    if ch
        gms = gmprint(A,B);
        str = [str sprintf('&\\begin{gmatrix}[p]\n') gms];
    end
    if j > 1
        if ch
            str = [str sprintf('\\rowops')];
        end
        ch = false;
        for k=j-1:-1:1
            multip = -A(k,j)/A(j,j);
            
            if multip ~= 0
                ch = true;
            else
                continue
            end
            
            A(k,:) = multip*A(j,:) + A(k,:);
            B(k,:) = multip*B(j,:) + B(k,:);
            
            str = [str sprintf('\n\\mult{%d}{%s}\n\\add{%d}{%d}', j-1,['\cdot ' sprintf(format_mult(multip))],j-1,k-1)];
        end
        if ch
            str = [str sprintf('\n\\end{gmatrix}\n\\\\ \n \\implies \n')];
        end
    end
end

% Normalize
if ch
    str = [str sprintf('\\rowops')];
end
ch = false;
for k = 1:r
    if A(k,k) ~= 1
        ch = true;
        nrm = 1/A(k,k);
        A(k,:) = A(k,:)*nrm;
        B(k,:) = B(k,:)*nrm;
        str = [str sprintf('\n\\mult{%d}{%s}', k-1,['\cdot ' sprintf(format_mult(nrm))])];       
    end
end
if ch
    str = [str sprintf('\n\\end{gmatrix}\n\\\\ \n\\implies\n')];
    str = [str sprintf('&\\begin{gmatrix}[p]\n') gmprint(A,B) '\end{gmatrix}'];
end


clipboard('copy', str)
disp('LaTeX code copied to clipboard.')
end

function ele_str = format_element(ele)
neg = '';
c = char(ele);
if c(1)=='-'
    neg = '-';
    ele = -ele;
end

str = sprintf('%s',ele);
parts=regexp(str,'([^/]*)/([^/]*)','tokens','once');
if length(parts) < 1
    ele_str = [neg str];
elseif parts{2}==1
    ele_str = [neg parts{1}];
else
   if strcmp(parts{2}([1,end]),'()')
       parts{2} = parts{2}(2:end-1);
   end
    ele_str = [neg '\\frac{' num2str(parts{1}) '}{' num2str(parts{2}) '}'];
end
end

function str = gmprint(A,B)
str = '';
for i = 1:size(A,1)
    for j = 1:size(A,2)
        str = [str format_element(A(i,j)) ' & '];
    end
    str = [str '| & '];
    for j = 1:size(B,2)
        str = [str format_element(B(i,j)) ' & '];
    end
    str = [str(1:end-2) '\\\\\n'];
end
str = sprintf([str(1:end-6) '\n']);
end

function ele_str = format_mult(ele)
neg = 0;
c = char(ele);
if c(1)=='-'
    neg = 1;
end

ele_str = format_element(ele);

if neg
    ele_str = ['\\left( ' ele_str ' \\right)'];
end
end