a=-3;
b=3;
for i=1:50
    r = a + (b-a).*rand(20,1);
    hist(r)
end