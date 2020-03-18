function r=Any_Radix_DIT(x,radix)
    N=length(x);
    M=ceil(N/radix);
    x1=zeros(1,radix*M);%Creating a zero squence of required length.
    x1(1:length(x))=x;%Superimposing the given sequence, so that the output is zero-padded.
    %Padding the Sequence.
   %%%%%%%%%%%%%%%%%%%%%%%%%Verified.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if (N>=radix)
       xmat=zeros(radix,M);%A blank matrix.
       for jr=1:M%This loop arranges 'x', column-wise.
           for ir=1:radix
               xmat(ir,jr)=x1(ir+(jr-1)*radix);
           end
       end
       xmat
       %Arranged Column-wise.
       %%%%%%%%%%%%%%%%%%%%%%%%%Verified.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       temp=xmat(1,:);
       step1=Any_Radix_DIT(temp,radix);
       for k=2:radix
           temp=xmat(k,1);
           for t=2:M
               temp=[temp,xmat(k,t)];
           end
           step1=[step1;Any_Radix_DIT(temp,radix)];
       end
       step1
       %Row-Wise DFT Achieved.
       %%%%%%%%%%%%%%%%%%%%%%%%Verified.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       phase=zeros(radix,M);
       for ir=1:radix
           for jr=1:M
               phase(ir,jr)=(ir-1)*(jr-1);
           end
       end
       step2=step1.*exp(-1i*phase*2*pi/N)
       %Phase matrix Multiplied.
       %%%%%%%%%%%%%%%%%%%%%%%%Verified.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       temp=reshape(step2(:,1),[1,radix]);
       step3=(fft(temp));
       for k=2:M
           temp=reshape(step2(:,k),[1,radix]);
           step3=[step3;fft(temp)];
       end
       step3=conj(step3)'
       %Column-Wise DFT Achieved.
       %%%%%%%%%%%%%%%%%%%%%%%%Verified.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       step3t=step3';
       r=step3t(:)';
   else
       r=fft(x)
   end
end
