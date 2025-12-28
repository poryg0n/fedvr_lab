     module util
     implicit none 
     contains

      subroutine open_files_to_store(workdir)
      implicit none
      include 'open/files_index'
      integer :: plot

       character(255) :: cwd
       character(3) :: str_run_id
       character(255) :: data_dir="data/"
       character(255) :: workdir
       character(128) :: command1
       character(128) :: command2
 

       include 'common/open_files'
       include 'open/files_specs'


       workdir="run_"

       if (plot.eq.0) then
          command2 = "ls -l data/ | grep -c ^d > data/run_id.txt"

          call getcwd(cwd)
          write(*,*) trim(cwd)

          call execute_command_line(command2, &
                    exitstat=status)

          open(1000,file='data/run_id.txt',status='unknown')
          read(1000,*) run_id
          run_id = run_id + 1

          write(str_run_id, "(i2.2)") run_id
          workdir=trim(data_dir)//trim(workdir)
          workdir=trim(workdir)//trim(str_run_id)//"/"
          print *, trim(workdir)
          command1 = "mkdir "//trim(workdir)
          call execute_command_line(command1, &
                    exitstat=status)
          call execute_command_line(command2, &
                    exitstat=status)
       end if

       include 'open/files'


       return
       end subroutine


       subroutine open_files_to_plot(arg)
       implicit none
       include 'open/files_index'
 
       character(3) :: arg
       character(255) :: data_dir="data/"
       character(255) :: workdir="run_"

 
       include 'common/open_files'
       include 'open/files_specs'
 
       workdir=trim(data_dir)//trim(workdir)
       workdir=trim(workdir)//trim(arg)//"/"
       print *, trim(workdir)

 
! ***  open stored data file for reading
       include 'open/files'
! ***  open projection data file for writing
       include 'open/files_processed'

  
       return
       end subroutine

      end
