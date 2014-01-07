%% OriTestEnv - created for MTex 3.3.2
% this script can be utilized to check consistency of MicroGenV1/COReV3
% mymath orientation functions with MTex

%% clear workspace
clc;
clear;

%% define filename

testid = 5;
%if (testid == 0) randomOrientation - test validity with import_core6orifile_v04.m
if (testid == 1) filename = 'final1.csv'; end %eul2quat
if (testid == 2) filename = 'final2.csv'; end %euleul2quat
if (testid == 3) filename = 'final3.csv'; end %rotateBaby
if (testid == 4) filename = 'final4.csv'; end %iamAquatDistance
if (testid == 5) filename = 'final5.csv'; end %eulscatter
    

%% import output data csv vom oritestenv data

newData1 = importdata(filename);

% Create new variables in the base workspace from those fields.
vars = fieldnames(newData1);
for i = 1:length(vars)
    assignin('base', vars{i}, newData1.(vars{i}));
end

clearvars i newData1 textdata vars

%% perform a number of tests on these data

%% test eul2quat2eul_luis (long n)
%//generate n random orientations, convert them into their quaternions and in Euler angle again
%//calculate occuring misorientation during each conversion step with misorientationCubic, ..quaternionCubic, 
%//cout<<"ID;ID;phi1;PHI;phi2;q0;q1;q2;q3;phi1;PHI;phi2;q0B;q1B;q2B;q3B;misorientationCubicEulerAB;misoriQxQAB;misoriQuatCubic0;misoriQuatCubic1;misoriQuatCubic2;misoriQuatCubic3" << endl;

if (testid == 1) 
    'eul2quat2eul'
    
    n = length(data(:,1)); 
    e = length(data(1,:));
    compdata = zeros(n,e);

    for i=1:n   
        qA = euler2quat(data(i,2), data(i,3), data(i,4),'Bunge');
        compdata(i,5) = get(qA,'a');
        compdata(i,6) = get(qA,'b');
        compdata(i,7) = get(qA,'c');
        compdata(i,8) = get(qA,'d');
        
        euler = Euler(qA,'Bunge');
        compdata(i,9) = euler(1);
        compdata(i,10) = euler(2);
        compdata(i,11) = euler(3);
        
        qB = euler2quat(compdata(i,9),compdata(i,10),compdata(i,11),'Bunge');
        compdata(i,12) = get(qB,'a');
        compdata(i,13) = get(qB,'b');
        compdata(i,14) = get(qB,'c');
        compdata(i,15) = get(qB,'d');     
        
        ori_qA = orientation('quaternion',compdata(i,5),compdata(i,6),compdata(i,7),compdata(i,8),symmetry('m-3m'));
        ori_qB = orientation('quaternion',compdata(i,12),compdata(i,13),compdata(i,14),compdata(i,15),symmetry('m-3m'));
        
        compdata(i,16) = angle(ori_qA, ori_qB);   
        
        if mod(i,100) == 0 
            i 
        end
    end    
end

save('final1.mtex.txt','compdata','-ascii');
csvwrite('final1.mtex.csv',compdata);
%% test eulereuler2disori_luis (long n)
%//generate n pairs of two random Euler angles, calculate disorientation between them
%//cout << "ID;ID;phi1A;PHIA;phi2A;phi1B;PHIB;phi2B;q0A;q1A;q2A;q3A;q0B;q1B;q2B;q3B;misorientationCubicEulerAB;misoriQxQAB;;misoriQuatCubic0;misoriQuatCubic1;misoriQuatCubic2;misoriQuatCubic3" << endl;

if (testid == 2) 
    'eulereuler2disori'
    
    n = length(data(:,1)); 
    e = length(data(1,:));
    disori_euler = zeros(n,1);
   
    for i=1:n 
        eA = Euler(data(i,2), data(i,3), data(i,4),'Bunge');
        eB = Euler(data(i,5), data(i,6), data(i,7),'Bunge');
        
        ori_eA = orientation(eA,symmetry('m-3m'));
        ori_eB = orientation(eB, symmetry('m-3m'));
                
        disori_euler(i) = angle(ori_eA,ori_eB); %yields same result: angle_outer(eA,eB);
        
        if mod(i,100) == 0 
            i 
        end
    end    
end

save('final2.euler.mtex.txt','disori_euler','-ascii');
csvwrite('final2.euler.mtex.csv',disori_euler);

%% test rotateBaby (long n) {
%//generates random Euler angle rotates that about random u,v,w axis about random angle
%//cout << "ID;ID;phi1A;PHIA;phi2A;u;v;w;randomAngle;eulRot[0];eulRot[1];eulRot[2]" << endl;

if (testid == 3) 
    'rotateBaby'
    
    n = length(data(:,1)); 
    e = length(data(1,:));
    
    compdata = zeros(n,3);
    for i=1:n 
        uvw = vector3d(data(i,5),data(i,6),data(i,7));
        angle = data(i,8);
        
        rotor = rotation('axis',uvw,'angle',angle);
           
        qA = euler2quat(data(i,2), data(i,3), data(i,4),'Bunge');
        qrotated = rotor * qA;
        [p1 P p2] = Euler(qrotated);
        
        compdata(i,1) = p1;
        compdata(i,2) = P;
        compdata(i,3) = p2;   
               
        if mod(i,100) == 0 
            i 
        end
    end    
end

save('final3.mtex.txt','compdata','-ascii');
csvwrite('final3.mtex.csv',compdata);

%% test iamAquatDistance long(n) {
%//generates two random euler angles, converts into quaternions and
%calculates distanceAB
%//cout << "ID;ID;phi1A;PHIA;phi2A;phi1B;PHIB;phi2B;q0A;q1A;q2A;q3A;q0B;q1B;q2B;q3B;distanceAB" << endl;

if (testid == 4) 
    'iamAquatAngle'
    
    n = length(data(:,1)); 
    e = length(data(1,:));
    compdata = zeros(n,e);

    for i=1:n   
        qA = euler2quat(data(i,2), data(i,3), data(i,4),'Bunge');
        compdata(i,8) = get(qA,'a');
        compdata(i,9) = get(qA,'b');
        compdata(i,10) = get(qA,'c');
        compdata(i,11) = get(qA,'d');
         
        qB = euler2quat(data(i,5),data(i,6),data(i,7),'Bunge');
        compdata(i,12) = get(qB,'a');
        compdata(i,13) = get(qB,'b');
        compdata(i,14) = get(qB,'c');
        compdata(i,15) = get(qB,'d');     
        
        %ori_qA = orientation('quaternion',compdata(i,8),compdata(i,9),compdata(i,10),compdata(i,11),symmetry('m-3m'));
        %ori_qB = orientation('quaternion',compdata(i,12),compdata(i,13),compdata(i,14),compdata(i,15),symmetry('m-3m'));
         
        compdata(i,16) = angle(qA, qB);       
        
        if (mod(i,100) == 0)
            i
        end
    end    
end

save('final4.mtex.txt','compdata','-ascii');
csvwrite('final4.mtex.csv',compdata);


%% test eulscatter

%//generates random euler angle and random scatter around that to form new orientation
%//cout<<"ID;ID;phi1A;PHIA;phi2A;myDev;eulNew0;eulNew1;eulNew2;disoriANew" << endl;
	
if (testid == 5) 
    'eulscatter'
    
    n = length(data(:,1)); 
    e = length(data(1,:));
    
    disori_anew = zeros(n,1);
   
    for i=1:n      
        qA = euler2quat(data(i,2), data(i,3), data(i,4),'Bunge');
        qNew = euler2quat(data(i,7), data(i,8), data(i,9),'Bunge');
           
        ori_qA = orientation(qA,symmetry('m-3m'));
        ori_qNew = orientation(qNew,symmetry('m-3m'));
        
        disori_anew(i) = angle(ori_qA, ori_qNew);
        
        if (mod(i, 100) == 0) 
            i
        end
    end    
end

save('final5.mtex.txt','disori_anew','-ascii');
csvwrite('final5.mtex.csv',disori_anew);



%% currently not tested!
%% test quat2rodrigues (long n) {
%//generates random euler angles, converts into quaternion and calculates Rodrigues type angle axis representation
%//cout<<"ID;ID;phi1A;PHIA;phi2A;q0A;q1A;q2A;q3A;angle;qr1;qr2;qr3" << endl;

if (testid == 4) 
    'quat2rodrigues'
    
    n = length(data(:,1)); 
    e = length(data(1,:));
    
    compdata = zeros(n,e);
    for i=1:n      
        qA = euler2quat(data(i,2), data(i,3), data(i,4),'Bunge');
        compdata(i,5) = get(qA,'a');
        compdata(i,6) = get(qA,'b');
        compdata(i,7) = get(qA,'c');
        compdata(i,8) = get(qA,'d');
        
        rodrig = Rodrigues(qA);
        
    
        rotor = rotation('axis',uvw,'angle',angle);
           
        qA = euler2quat(data(i,2), data(i,3), data(i,4),'Bunge');
        qrotated = rotor * qA;
        [p1 P p2] = Euler(qrotated);
        
        compdata(i,1) = p1 / 180 * pi;
        compdata(i,2) = P / 180 * pi;
        compdata(i,3) = p2 / 180 * pi;              
    end    
end

