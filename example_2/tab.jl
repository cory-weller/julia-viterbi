#5 major command-line inputs are required
#1: Test Chromosome
#2: Observation chromosome
#3: Text file used to obtain the appropriate columns
#4: File to obtain the file's size
#5: Specifies the number of founders that you want

using Plots
global start=0
global finish=100000
dataset=ARGS[1]                                              #Test chromosome
observations=ARGS[2]                                         #Observations
little=ARGS[4]                                              #File to obtain file size
size1=readlines(little)
size2=split.(size1,"\t")
global l=parse(Int,size2[1][2])
global chrom=size2[1][1]                                   #Chromosome that you are reading,  might need this

tabix1=`tabix -p vcf $dataset $chrom:$start-$finish`        #This is outside of the loop to fill the first row of total[]
tabix2=`tabix -p vcf $observations $chrom:$start-$finish`
x=readlines(tabix1)
x2=readlines(tabix2)
y=split.(x,"\t")
first=y[1]
global founders=length(first)-4
y2=split.(x2,"\t")
a=size(y)
b=size(y2)
global len=a[1]
global farray=zeros(founders)
for i in 1:1:len
    for j in 1:1:founders
        if y[i][j+4]==y2[i][3]
            farray[j]=farray[j]+1
        end
    end
end
global total=[farray[1]/len]
for i in 2:1:founders
    n=[farray[i]/len]
    global total=hcat(total,n)
end
start+=50000
finish+=50000

while finish<l
    tabix1=`tabix -p vcf $dataset $chrom:$start-$finish`
    tabix2=`tabix -p vcf $observations $chrom:$start-$finish`
    x=readlines(tabix1)
    x2=readlines(tabix2)
    y=split.(x,"\t")
    y2=split.(x2,"\t")
    a=size(y)
    b=size(y2)
    global len=a[1]
    global farray=zeros(founders)
    for i in 1:1:len
        for j in 1:1:founders
            if y[i][j+4]==y2[i][3]
                farray[j]=farray[j]+1
            end
        end
    end
    global con=[farray[1]/len]
    for i in 2:1:founders
        n=[farray[i]/len]
        con=hcat(con,n)
    end
    global total=vcat(total,con)
    global start+=50000
    global finish+=50000
end

#global cutoff=parse(Int,ARGS[5])
#cutoff=cutoff/10

siz=size(total)
global carray=zeros(founders)
for i in 1:1:siz[1]
    for j in 1:1:founders
        if total[i,j] >= 0.9                 #Cutoff for high concordance is set at 0.9
            carray[j]=carray[j]+1
        end
    end
end


global numfounds=parse(Int,ARGS[5])                 #How many founders you want to run the Viterbi algorithm on
global max=[]
for i in 1:1:numfounds
    global maxval=0
    global counter=0
    for j in 1:1:founders
        if carray[j]>maxval
            maxval=carray[j]
            counter=j
        end
    end
    push!(max,counter)
    carray[counter]=0
end

txtfile=ARGS[3]                                             #Actual text file used to obtain the appropriate columns

if numfounds==1
    m11=max[1]+4                                            #Add four to get the actual founder column
    text=`cut -f $m11 $txtfile`
end
if numfounds==2
    m11=max[1]+4                                            #Add four to get the actual founder column
    m12=max[2]+4                                            #Because of the four extra columns
    text=`cut -f $m11,$m12 $txtfile`
end
if numfounds==3
    m11=max[1]+4                                            #Add four to get the actual founder column
    m12=max[2]+4                                            #Because of the four extra columns
    m13=max[3]+4
    text=`cut -f $m11,$m12,$m13 $txtfile`
end
if numfounds==4
    m11=max[1]+4                                            #Add four to get the actual founder column
    m12=max[2]+4                                            #Because of the four extra columns
    m13=max[3]+4
    m14=max[4]+4
    text=`cut -f $m11,$m12,$m13,$m14 $txtfile`
end
if numfounds==5
    m11=max[1]+4                                            #Add four to get the actual founder column
    m12=max[2]+4                                            #Because of the four extra columns
    m13=max[3]+4
    m14=max[4]+4
    m15=max[5]+4
    text=`cut -f $m11,$m12,$m13,$m14,$m15 $txtfile`
end
if numfounds==6
    m11=max[1]+4                                            #Add four to get the actual founder column
    m12=max[2]+4                                            #Because of the four extra columns
    m13=max[3]+4
    m14=max[4]+4
    m15=max[5]+4
    m16=max[6]+4
    text=`cut -f $m11,$m12,$m13,$m14,$m15,$m16 $txtfile`
end
if numfounds==7
    m11=max[1]+4                                            #Add four to get the actual founder column
    m12=max[2]+4                                            #Because of the four extra columns
    m13=max[3]+4
    m14=max[4]+4
    m15=max[5]+4
    m16=max[6]+4
    m17=max[7]+4
    text=`cut -f $m11,$m12,$m13,$m14,$m15,$m16,$m17 $txtfile`
end
if numfounds==8
    m11=max[1]+4                                            #Add four to get the actual founder column
    m12=max[2]+4                                            #Because of the four extra columns
    m13=max[3]+4
    m14=max[4]+4
    m15=max[5]+4
    m16=max[6]+4
    m17=max[7]+4
    m18=max[8]+4
    text=`cut -f $m11,$m12,$m13,$m14,$m15,$m16,$m17,$m18 $txtfile`
end
run(text)

plotly(size=(600,400))
heatmap(total, xlabel="Founder", ylabel="Genomic Window", show=true, fmt= :png)
