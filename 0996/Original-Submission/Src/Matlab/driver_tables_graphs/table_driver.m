%
% Input solvePOParrow.csv and solvePOPchordal.csv
% 
% Output POPtable.txt --- The original data for POPtable1.txt and POPtable1.txt for 
%                         Tables 1 and 2 in Section 4.4 of TOMS's paper
%

percent = char(37); 
id = @(x) x;
fileID = fopen('POPdata.txt', 'w');
%%
fprintf(fileID, '\n\n\\def\\CC{\\mbox{%s\\cal C%s}}\n\n',char(36),char(36));
fprintf(fileID, '\n\n%s Arrow\n\n',percent);
ta = readtable('solvePOParrow.csv');
ta.box = double(~ta.bin);
ce = table2cell(ta);
%[box a b c d n l], [comp solver], [LBv time iter bpiter term]
out= pivottable(ce, [16 3 4 5 1 2 6], [8 10], [11 12 13 15 14], {id,id,id,id,id});
values = out(3:end, 5:end);
%%
offset = 8;
numsolvers = size(out(1, offset:end), 2) / 5;
fprintf(fileID, '\\begin{landscape}\n');
fprintf(fileID, '\\begin{table}\n');
fprintf(fileID, '\\scalebox{0.6}{\n');
fprintf(fileID, ['\\begin{tabular}{|c|r|r|' repmat('r@{~(}r@{,~}r@{:}r@{,~}r@{)~~~}|', 1, numsolvers) '}\\hline\n']);
fprintf(fileID, '  \\multicolumn{3}{|c|}{obj. $\\backslash$ constr.}');
for solblk = 0:(numsolvers-1)
    if ~out{1, offset + solblk * 5} % i.e., comp=0
        comp = '$\emptyset$';
    else
        comp = 'random';
    end
    fprintf(fileID, ' & \\multicolumn{5}{c|}{$\\CC:$ %s}', comp);
end
fprintf(fileID, ' \\\\ \\hline\n');
%
fprintf(fileID, ' & & ');
for solblk = 0:(numsolvers-1)
    solver = out{2, offset + solblk * 5};
    fprintf(fileID, ' & \\multicolumn{5}{c|}{%s}', solver);
end
fprintf(fileID, ' \\\\\n');
%
fprintf(fileID, 'type & d & n');
for solblk = 0:(numsolvers-1)
    solver = out{2, offset + solblk * 5};
    fprintf(fileID, ' & \\multicolumn{5}{r|}{LBv (sec,apgit:bpit,term)}');
end
fprintf(fileID, '\\\\ \\hline\n');
%%
[blkid, ~, ic] = unique(cell2mat(out(3:end, 1:4)), 'rows');
for blk = 1:size(blkid, 1)
    if ~blkid(blk, 1) % i.e., box=0 (see line 6.)
        v = 'bin';
    else
        v = 'box';
    end

    blk_mask = ic == blk;
    format = ['\\multirow{' num2str(sum(blk_mask)) '}{*}{\\shortstack{Arrow ' v '\\\\($a=%d$, \\\\$b=%d$, \\\\$c=%d$)}}\n'];
    fprintf(fileID, format, blkid(blk, 2:end));

    formatSpecIdx = ' & %d & %d($\\ell$=%d)';
    formatSpecVal = ' & %1.6e & %1.2e & %d & %d & %d';
    formatSpec = [formatSpecIdx repmat(formatSpecVal, 1, 4) '\\\\\n'];
    blkvalues = values(blk_mask, :);
    for row = 1:sum(blk_mask)
        fprintf(fileID, formatSpec, blkvalues{row, :});
    end
    fprintf(fileID, '\\hline\\hline\n');
end
fprintf(fileID, '\\end{tabular}\n');
fprintf(fileID, '}\n');
fprintf(fileID, '\\end{table}\n');
fprintf(fileID, '\\end{landscape}\n');

%%
fprintf(fileID, '\n\n%s Chordal \n\n',percent);
ta = readtable('solvePOPchordal.csv');
ta.box = double(~ta.bin);
ce = table2cell(ta);
%[box r d n], [comp solver], [LBv time iter bpiter term]
out = pivottable(ce, [13 3 1 2], [5 7], [8 9 10 12 11], {id,id,id,id,id});

values = out(3:end, 3:end);
[blkid, ~, ic] = unique(cell2mat(out(3:end, 1:2)), 'rows');
%%
offset = 5;
numsolvers = size(out(1, offset:end), 2) / 5;
fprintf(fileID, '\\begin{landscape}\n');
fprintf(fileID, '\\begin{table}\n');
fprintf(fileID, '\\scalebox{0.6}{\n');
fprintf(fileID, ['\\begin{tabular}{|c|r|r|' repmat('r@{~(}r@{,~}r@{:}r@{,~}r@{)~~~}|', 1, numsolvers) '}\\hline\n']);
fprintf(fileID, '  \\multicolumn{3}{|c|}{obj. $\\backslash$ constr.}');
for solblk = 0:(numsolvers-1)
    if ~out{1, offset + solblk * 5}  % i.e., comp=0
        comp = '$\emptyset$';
    else
        comp = 'random';
    end
    fprintf(fileID, ' & \\multicolumn{5}{c|}{$\\CC:$ %s}', comp);
end
fprintf(fileID, ' \\\\ \\hline\n');
%
fprintf(fileID, ' & & ');
for solblk = 0:(numsolvers-1)
    solver = out{2, offset + solblk * 5};
    fprintf(fileID, ' & \\multicolumn{5}{c|}{%s}', solver);
end
fprintf(fileID, ' \\\\\n');
%
fprintf(fileID, 'type & d & n');
for solblk = 0:(numsolvers-1)
    solver = out{2, offset + solblk * 5};
    fprintf(fileID, ' & \\multicolumn{5}{r|}{LBv (sec,apgit:bpit,term)}');
end
fprintf(fileID, '\\\\ \\hline\n');
%%
for blk = 1:size(blkid, 1)
    if ~blkid(blk, 1)  % i.e., box=0
        v = 'bin';
    else
        v = 'box';
    end

    blk_mask = ic == blk;
    format = ['\\multirow{' num2str(sum(blk_mask)) '}{*}{\\shortstack{Chordal ' v '\\\\($r=%1.1f$)}}\n'];
    fprintf(fileID, format, blkid(blk, 2:end));

    formatSpecIdx = ' & %d & %d';
    formatSpecVal = ' & %1.6e & %1.2e & %d & %d & %d';
    formatSpec = [formatSpecIdx repmat(formatSpecVal, 1, 4) '\\\\\n'];
    blkvalues = values(blk_mask, :);
    for row = 1:sum(blk_mask)
        fprintf(fileID, formatSpec, blkvalues{row, :});
    end
    fprintf(fileID, '\\hline\\hline\n');
end
fprintf(fileID, '\\end{tabular}\n');
fprintf(fileID, '}\n');
fprintf(fileID, '\\end{table}\n');
fprintf(fileID, '\\end{landscape}\n');
%%
fclose(fileID);