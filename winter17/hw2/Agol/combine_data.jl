# Combine Keck, SOPHIE, ELODIE datasets (1=ELODIE; 2=HIRES; 3=SOPHIE)

data1 = readdlm("mystery_planet1a.txt")
data2 = readdlm("mystery_planet1b.txt")
data3 = readdlm("mystery_planet1c.txt")
data3[:,2] /= 1000.
data3[:,3] /= 1000.

period = 111.43670
clf()
offset1 = [3.767*ones(61)]
scatter(mod(data1[:,1],period),data1[:,2]-offset1)
offset2 = [3.788*ones(19);3.911*ones(48)]
scatter(mod(data2[:,1],period),data2[:,2]-offset2,color="Red")
scatter(mod(data3[:,1],period),data3[:,2],color="Green")

n1 = 61
n2 = 67
n3 = 46
datatot=zeros(n1+n2+n3,3)
for i=1:n1
  datatot[i,1]=data1[i,1]
  datatot[i,2]=data1[i,2] -offset1[i]
  datatot[i,3]=data1[i,3]
end
for i=1:n2
  datatot[i+n1,1]=data2[i,1]
  datatot[i+n1,2]=data2[i,2] -offset2[i]
  datatot[i+n1,3]=data2[i,3]
end
for i=1:n3
  datatot[i+n1+n2,1]=data3[i,1]
  datatot[i+n1+n2,2]=data3[i,2]
  datatot[i+n1+n2,3]=data3[i,3]
end


writedlm("mystery_planet1.txt",datatot)
