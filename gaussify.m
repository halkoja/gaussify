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
nl = sprintf('\n');

[r,~] = size(A);

ch = true;
% Downward Elimination
for j = 1:r
    if ch
        gms = gmprint(A,B);
        str = [str '&\begin{gmatrix}[p]' nl gms];
    end
    if j < r
        if ch
            str = [str '\rowops' nl];
        end
        ch = false;
        for k=j+1:r
            multip = -A(k,j)/A(j,j);
            
            if multip ~= 0
               ch = true;
            else
                continue
            end
            
            A(k,:) = simplifyFraction(multip*A(j,:) + A(k,:));
            B(k,:) = simplifyFraction(multip*B(j,:) + B(k,:));
            
            str = [str sprintf('\\mult{%d}{%s}\n\\add{%d}{%d}\n', j-1,['\cdot ' format_mult(multip)],j-1,k-1)];
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
            
            A(k,:) = simplifyFraction(multip*A(j,:) + A(k,:));
            B(k,:) = simplifyFraction(multip*B(j,:) + B(k,:));
            
            str = [str ...
                sprintf('\\mult{%d}{%s}\n\\add{%d}{%d}\n', ...
                j-1,['\cdot ' format_mult(multip)],j-1,k-1)];
        end
        str = [str '\end{gmatrix}' nl];
        if ch
            str = [str '\\' nl '\implies'];
        end
    end
end

% Upwards
for j = r-1:-1:1
    if ch
        gms = gmprint(A,B);
        str = [str '&\begin{gmatrix}[p]' nl gms];
    end
    if j > 1
        if ch
            str = [str '\rowops'];
        end
        ch = false;
        for k=j-1:-1:1
            multip = -A(k,j)/A(j,j);
            
            if multip ~= 0
                ch = true;
            else
                continue
            end
            
            A(k,:) = simplifyFraction(multip*A(j,:) + A(k,:));
            B(k,:) = simplifyFraction(multip*B(j,:) + B(k,:));
            
            str = [str sprintf('\n\\mult{%d}{%s}\n\\add{%d}{%d}', j-1,['\cdot ' format_mult(multip)],j-1,k-1)];
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
        A(k,:) = simplifyFraction(A(k,:)*nrm);
        B(k,:) = simplifyFraction(B(k,:)*nrm);
        str = [str sprintf('\n\\mult{%d}{%s}', k-1,['\cdot ' format_mult(nrm)])];       
    end
end
if ch
    str = [str sprintf('\n\\end{gmatrix}\n\\\\ \n\\implies\n')];
    str = [str sprintf('&\\begin{gmatrix}[p]\n') gmprint(A,B) '\end{gmatrix}'];
end


clipboard('copy', str)
disp('LaTeX code copied to clipboard.')
end

function str = gmprint(A,B)
str = '';
for i = 1:size(A,1)
    for j = 1:size(A,2)
        str = [str latex(A(i,j)) ' & '];
    end
    str = [str '| & '];
    for j = 1:size(B,2)
        str = [str latex(B(i,j)) ' & '];
    end
    str = [str(1:end-2) '\\' sprintf('\n')];
end
str = [str(1:end-4) sprintf('\n')];
end

function ele_str = format_mult(ele)
ele_str = latex(ele);

if ele_str(1) == '-'
    ele_str = ['\left( ' ele_str ' \right)'];
end
end