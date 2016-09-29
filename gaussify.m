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
gbeg = ['&\begin{gmatrix}[p]' nl];
gend = '\end{gmatrix}';
imp = ['\implies' nl];

[r,~] = size(A);

ch = true;
% Downward Elimination
for j = 1:r-1
    if ch
        gms = gmprint(A,B);
        str = [str gbeg gms];
    end
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
        
        str = [str mulstr(j-1,multip) addstr(j-1,k-1)];
    end
    if ch
        str = [str gend '\\' nl imp];
    end
end

% Upwards
for j = r:-1:2
    if ch
        gms = gmprint(A,B);
        str = [str gbeg gms];
    end
    if j > 1
        if ch
            str = [str '\rowops' nl];
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
            
            str = [str mulstr(j-1,multip) addstr(j-1,k-1)];
        end
        if ch
            str = [str gend '\\' nl imp];
        end
    end
end

% Normalize
if ch
    str = [str gbeg gmprint(A,B) '\rowops' nl];
end
ch = false;
for k = 1:r
    if A(k,k) ~= 1
        ch = true;
        nrm = 1/A(k,k);
        A(k,:) = simplifyFraction(A(k,:)*nrm);
        B(k,:) = simplifyFraction(B(k,:)*nrm);
        str = [str mulstr(k-1,nrm)];
    end
end
if ch
    str = [str gend '\\' nl imp];
    str = [str gbeg gmprint(A,B) gend];
else
    str = [str gend];
end


clipboard('copy', str)
disp('LaTeX code copied to clipboard.')
end

function str = gmprint(A,B)
splA = tidytex(A);
splB = tidytex(B);

spl = [splA splB];
spl = arrayfun(@(x) strjoin(spl(x,:),' &|& '), 1:size(spl,1),'UniformOutput',false);
str = strjoin(spl,' \\\\\n');
str(end+1) = sprintf('\n');
end

function spl = tidytex(A)
strA = latex(A);
strA = strrep(strA,'\end{array}\right)','');
strA = regexprep(strA,'\\left(\\begin{.*}{c+}','');
spl = strsplit(strA,'\\\\')';
end

function ele_str = format_mult(ele)
ele_str = latex(ele);

if ele_str(1) == '-'
    ele_str = ['\left( ' ele_str ' \right)'];
end
end

function astr = addstr(r1,r2)
    astr = ['\add{' num2str(r1) '}{' num2str(r2) '}' sprintf('\n')];
end
function mstr = mulstr(r,mult)
    mstr = ['\mult{' num2str(r) '}{\cdot ' num2str(format_mult(mult)) '}' sprintf('\n')];
end