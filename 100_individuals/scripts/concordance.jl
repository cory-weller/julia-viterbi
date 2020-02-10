#6 major command-line inputs are required
#1: Test Chromosome
#2: Observation chromosome
#3: Text file used to obtain the appropriate columns
#4: File to obtain the file's size
#5: Specifies the number of founders that you want
#6: Cutoff value

####Look at comment telling you next step

#using Plots
global start=0
global finish=200000
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
#global founders=length(first)-4
global founders=length(first)-9                          ######!!!!!!!!!This will change depending on the VCF file
y2=split.(x2,"\t")
a=size(y)
b=size(y2)
global kant1=[]
global kant1count=0
global len=b[1]
global farray=zeros(founders)
for j in 1:1:founders
    global kant1count=0
    for i in 2:1:len
        if y2[i][3]!="."
            global kant1count+=1                   ###### Must make sure that y and y2 are the same size, so y's size must be changed
            if y[i][j+9]==y2[i][3]                         ######!!!!!This will change!!!
                farray[j]=farray[j]+1
            end
        end
    end
    push!(kant1,kant1count)
end
global total=[farray[1]/kant1[1]]
for i in 2:1:founders
    n=[farray[i]/kant1[i]]
    global total=hcat(total,n)
end
start+=100000
finish+=100000

while finish<l
    tabix1=`tabix -p vcf $dataset $chrom:$start-$finish`
    tabix2=`tabix -p vcf $observations $chrom:$start-$finish`
    x=readlines(tabix1)
    x2=readlines(tabix2)
    y=split.(x,"\t")
    y2=split.(x2,"\t")
    a=size(y)
    b=size(y2)
    global kant=[]
    global kantcount=0
    global len=b[1]
    global farray=zeros(founders)
    for j in 1:1:founders
        global kantcount=0
        for i in 2:1:len                                                  ####Starts at 2 because of unnescessary first row
            if y2[i][3]!="."
                global kantcount+=1
                if y[i][j+9]==y2[i][3]                         ######!!!!!This will change!!!
                    farray[j]=farray[j]+1
                end
            end
        end
        push!(kant,kantcount)
    end
    global con=[farray[1]/kant[1]]
    for i in 2:1:founders
        n=[farray[i]/kant[i]]
        con=hcat(con,n)
    end
    global total=vcat(total,con)
    global start+=50000
    global finish+=50000
end

global cutoff=parse(Int,ARGS[6])
if cutoff<10
    cutoff=cutoff/10
else
    cutoff=cutoff/100
end

siz=size(total)
global carray=zeros(founders)               ##carray is count array
for i in 1:1:siz[1]
    for j in 1:1:founders
        if total[i,j] >= cutoff                 #Cutoff for high concordance is set at 0.9
            carray[j]=carray[j]+1
        end
    end
end

#global numfounds=parse(Int,ARGS[5])                 #How many founders you want to run the Viterbi algorithm on
global numfounds=16
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
    if counter!=0
        carray[counter]=0
    end
end

tf=ARGS[5]
tf2=readlines(tf)
sizzle1=size(tf2)[1]
parsed=[]
for i in 1:1:sizzle1
    num=parse(Int,tf2[i])-9                   ####Must be subtracted so the rows match, subject to change
    push!(parsed,num)
end
sizzle2=size(max)[1]
global finds=0
for i in 1:1:sizzle1
    for j in 1:1:sizzle2
        if parsed[i]==max[j]
            global finds=finds+1
        end
    end
end
println(finds/sizzle1)
print(max)

txtfile=ARGS[3]                                             #Actual text file used to obtain the appropriate columns
newarray=split.(txtfile,"\t")

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
#run(text)

#plotly(size=(600,400))
#heatmap(total, xlabel="Founder", ylabel="Genomic Window", show=true, fmt= :png)
