% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.

% First derivative flux image analysis used in NCOMM 2021 (Mo)

function imroldiff(path,w)
% seQ: image data
% imread_big/imwrite_big: see code distributed by YoonOh Tak, GIST, Copyright (c) 2012
tic;
temp=strsplit(path,'\'); folderpath=path(1:strfind(path,temp{end})-1);
filename=strsplit(temp{end},'.'); filename=filename{1};
seQ=imread_big(path);
fprintf('%s%s\n',folderpath,filename); 
kL=zeros([size(seQ,1:2),size(seQ,3)-w+1]); mL=kL;
for i=1:size(seQ,1)
    for j=1:size(seQ,2)
    [kL(i,j,:),mL(i,j,:)]=rollingDiff(seQ(i,j,:),w);
    end
end
kL=single(kL); mL=single(mL);
xL=single(kL.*(mL/max(mL,[],'all')));
imwrite_big(kL,[folderpath filename '_diff_' num2str(w) '.tif']);
imwrite_big(mL,[folderpath filename '_mean_' num2str(w) '.tif']);
imwrite_big(xL,[folderpath filename '_dscale_' num2str(w) '.tif']);
toc;


function [d,m]=rollingDiff(seq,w)
m=[]; % Rolling mean
for i=1:(length(seq)-w+1)
    t=seq(i:(i+w-1));
    m=[m mean(t)];
end
d=[0 diff(m)];
