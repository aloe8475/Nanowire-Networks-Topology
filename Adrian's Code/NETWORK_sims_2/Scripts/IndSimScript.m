% script for running independent simmulations.

% Network structure: 

% Settings structure of Network structure. (Network.Settings).
% any parameter is iterable over different simulations with same network
% and electrodes just by overrriding that parameter in settings structure
% in the control loop of the script. 

  
    
        

        
       
        
        %% open file
        
        [filename, pathname] = uigetfile('*.mat', 'Select the network');
        if isequal(filename,0)
            disp('User selected Cancel')
            return;
        else
            disp(['User selected ', fullfile(pathname, filename)])
        end
        
        NetS=load(fullfile(pathname,filename),'-mat');
        
        CurrNet=NetS.Network;  

        
        %% All simulations settings (EXCLUDING ELECTRODES POSITION AND SIGNAL)
        sim=struct();
        sim.PreName='Sim';
        sim.ElectrodeType='point';
        sim.OneNW=1;
        sim.SquareNW=0;
        
        sim.IniRes=5e6;
        sim.Time=1;
        sim.Step=0.01;
        sim.Ron=5e3;
        sim.Roff=5e6;
        sim.MaxW=5e-9;
        sim.Mobility=0.5e-12;
        sim.Tau=1;
        sim.SigmaW=0.1;
        sim.SigmaNoise=0.1;
        sim.IniW=1e-11;
        sim.Alpha=1e-3;
        sim.Factor=sim.Mobility*sim.Ron/sim.MaxW;
        sim.MaxI=1e-4;
        sim.WForm=5e-9;
        sim.FactorZ=1;
        sim.WDissolve=1e-9;
        sim.Random=0;
        sim.Break=1;
        sim.ThreshBreak=1e-9;
        
        sim.Model='HP';
        sim.Pow=0;
        
        
        sim.ElectrodesInfo={};
        sim.Vmax=1;
        sim.Vmin=1e-3;
        sim.NoC=1;
        sim.Frequency=0;
        sim.Name='Source1';
        sim.OpenCheck=0;
        sim.Duty=50;
        sim.SigType='Constant';
        sim.TStartOff=0;
        sim.SetFreq=10;
        sim.VSet=5;
        sim.ReadFreq=1;
        sim.VRead=5;
        sim.VStartOff=0;
        sim.ReadDuty=5;
        sim.SetDuty=5;
        sim.SType=1;
        sim.DType=0;
        
        %% electrode settings
         %as electrodes and time are already set in the main gui, it is dangerous
        %to change here the position of the electrodes, but maybe the
        % type of waveform send to the electrodes is feasible, or the duty
        %cycle of the waveform, or duration.
        % You can change the basic electrode configuration parameters here
        %Similar as in the sim conf. GUI. And call the ChangeElectrodes
        %function afterward with the electrode(s) you want to change.
        % Electrodes is a table object contained within the 
        %Simulation struct. 
        
       
        
        
        elsettings.Time=CurrNet.Simulations{1}.Settings.Time;
        elsettings.Step=CurrNet.Simulations{1}.Settings.Step;
        elsettings.Vmax=0.5;
        elsettings.Vmin=0.001;
        elsettings.NoC=1;
        elsettings.Frequency=1;
        elsettings.Duty=0.5;
         %Types: {'Constant','Triangular(IV)','Square','SquarePWM'}
        elsettings.SigType='Constant';
        %only for squarepwm settings
        elsettings.TStartOff=0;
        elsettings.SetFreq=1;
        elsettings.ReadFreq=5;
        elsettings.VSet=5;
        elsettings.VRead=5;
        elsettings.VStartOff=0.1;
        elsettings.NoP=elsettings.Time/elsettings.Step;
        elsettings.OpenCheck=0;
        
        %Number of electrodes
        Nel=height(CurrNet.Simulations{1}.Electrodes);
        
        %Name of electrodes
        Names=CurrNet.Simulations{1}.Electrodes.Name;
        %change electrode number 1
        [CurrNet.Simulations{1}.Time,CurrNet.Simulations{1}.Electrodes]...
            =ChangeElectrodes(CurrNet.Simulations{1}.Electrodes,1,elsettings);
        
        %% preparing the loop
        
        %parameter you might want to change, as R off, or V in the
        %electrodes
        
        %Roff=linspace(1e6,10e6,3);
        Vmax=linspace(0.4,1,3);
       
        % NUMBER OF ITERATIONS
        n_iter=length(Roff);
        
        % Warning: Simulation parameters and electrode info copied form GUI
        %
        NewSim=CurrNet.Simulations{1};
        
        %updating with sims structure of this script (optional)
        
        %NewSim.Settings=sim;
        
        % add (1) or copy (2) simulation
        add_flag=1;
        
        % incase add_flag==2;
        if isequal(add_flag,2)
            Networks=cell(1,n_iter);
        end
        
    
        
        
        %%loop
        for i=1:n_iter
            
            %here you can change parameter you want to override
            % i put by default Roff,but almost everything in sims structure
            %, should be fine.
            
            %NewSim.Settings.Roff=Roff(i);
           
            %here's an example of changing a parameter for the electrodes,
            %once you change the parameter, you need to call to the
            %ChangeElectrodes function, in this case, i keep under comments
            %the code to change the constant voltage from 0.5 to 1.
           
            elsettings.Vmax=Vmax(i);
            [NewSim.Time,NewSim.Electrodes]...
            =ChangeElectrodes(NewSim.Electrodes,1,elsettings);
            
            
            
            %simulation
            tic;
            [NewSim.Data,NewSim.Name]=LaunchSimulation(NewSim);
            
            
            %updating strucutres and time
            disp(strcat('Simulation completed. Elapsed time is',{' '},num2str(toc),{' '},'seconds'));
            
            NewSim.SimInfo=GetSimInfo(NewSim.Data);
            
            %storing on cell array of simulations under CurrNet.Simulations
            %field or copying network directly
            if isequal(add_flag,1)
                CurrNet.Simulations{i}=NewSim;
            else
                Networks{i}=CurrNet;
                Networks{i}.Simulations{1}=NewSim;
            end
            
            
        end
        
        
        
        %% save data
        
        %add new name information as you wish, i just put SC for
        %scripted, but maybe you may want to add relevant information
        
        %you can later open these files on the main GUI using the open
        %network function
        
        
        % path for saved networks
        
        save_path='C:\Users\Nanofig\Desktop\ScriptNets\';
        
        if isequal(add_flag,1)        
            CurrNet.Name=strcat(CurrNet.Name,'_SC');
            Name=CurrNet.Name;
            Name=strrep(Name,':','_');
            Name=strrep(Name,'/','_');
            Name=strcat(Name,'.mat');
            
            %in main_sims, i use the struct named SelNet to load saved
            %networks, thats why i copy the struct here with that name
            SelNet=CurrNet;
            save(fullfile(save_path,Name),'SelNet');
        else
            for i=1:length(Networks)
                Networks{i}.Name=strcat(Networks{i}.Name,'_SC_',num2str(i));
                Name=Networks{i}.Name;
                Name=strrep(Name,':','_');
                Name=strrep(Name,'/','_');
                Name=strcat(Name,'.mat');
                SelNet=Networks{i};
                save(fullfile(save_path,Name),'SelNet');
            end          
        end
        
        
        %% clear all variables from workspace?
        
         % clear all
        
        
        
        
        
        
        
        
        
        