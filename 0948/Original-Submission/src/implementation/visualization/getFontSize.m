function fontSize = getFontSize(n)

if n <= 20, fontSize = 10;
elseif n <= 30, fontSize = 8;
elseif n <= 40, fontSize = 7;
elseif n <= 55, fontSize = 6;
elseif n <= 70, fontSize = 5;
else fontSize = 4;
end
end