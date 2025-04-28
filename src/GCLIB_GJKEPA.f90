    !  ͨ�ñ����Fortran���ݽṹ��������  (GCLIB)
    !   General and Convenient Fortran Data Structure and Function Library
    !---------------------------------------------------------------
    !
    !       ��GJK - EPA ��ײ����㷨��
    !       Gilbert-Johnson-Keerthi  -   Expanding Polytope Algorith
    !       ( ��������-Լ��ѷ-����ϣ�㷨 - ��չ�������㷨 )
    !       2023.5.9 - code by XIE Jihong
    !       2023.9.7 - omp optimization by XIE Jihong
    !       ʹ�÷����� ֱ�ӵ���GJKEPA()���ɻ����ײ��Ϣ
    !---------------------------------------------------------------
    MODULE GCLIB_GJKEPA
    USE GCLIB_List, ONLY : List_Array3d, ListNode_Array3d
    USE GCLIB_QuickHull, ONLY : QuickHull
    USE GCLIB_DeHull, ONLY : getHullMeshesVertex, getHullMeshesVertex_post
    USE OMP_LIB
    IMPLICIT NONE
        !---------------------------------------------------------------
        ! GJK 
        PUBLIC :: GJKEPA
        PRIVATE :: EPA_solu
        PRIVATE :: get_collisionPoint_01, get_collisionPoint_02
        PRIVATE :: get_info_collisionType, get_nearest_points

        !PUBLIC :: GJK
        PRIVATE :: support_mapping, update_simplex_GJK
        PRIVATE :: RoughCollisionDetection_SphericalEnvelope
        
        !---------------------------------------------------------------
        ! Random data Lib
        PRIVATE :: GET_RANDOM_UNIT_VECTOR
        
        !---------------------------------------------------------------
        ! MATH TOOLS
        PRIVATE :: CROSS_PRODUCT_3D, isPointInSimplex, IS_INSIDE_PF, UTZVEC
        PRIVATE :: DIST_PF_SIGN, UNINML, OVERLAP, VEC_PL, FOOT_LL, FOOT_PL, SORT_CLOCK
        
    CONTAINS
        SUBROUTINE GJKEPA(version_, TOL_FF_, &
            p1_, p2_, collision_, colliType_, &
            nearest_points_, collision_normal_, collision_point_, penetration_depth_)
        !DEC$ ATTRIBUTES DLLEXPORT :: GJKEPA
            IMPLICIT NONE
            INTEGER*4, INTENT(IN) :: version_ ! �汾��ǣ���ǰ����1 ��2����Ӱ�������õĺ�����
	        REAL*8, INTENT(IN) :: TOL_FF_   ! ͨ��=1�����ֵԽ���ж�Ϊ����Ӵ���������Խ���ɡ�
            REAL*8, INTENT(IN) :: p1_(:,:), p2_(:,:)  ! �������������͹������Ķ��㼯
            LOGICAL*1, INTENT(OUT) :: collision_ ! �Ƿ���ײ
            INTEGER*4, INTENT(OUT) :: colliType_ ! ��ײ���ͣ�0 - δ��ײ�� 1 - ������ײ���� �� 2 - ������ײ(��ֵ�ж�)
            REAL*8, INTENT(OUT) :: nearest_points_(2,3)     ! ������
            REAL*8, INTENT(OUT) :: collision_normal_(3)     ! ��ײ����
            REAL*8, INTENT(OUT) :: collision_point_(3)      ! ��ײ��
            REAL*8, INTENT(OUT) :: penetration_depth_       ! ��͸���
            !---------------------------------------------------------------
            REAL*8, SAVE :: simplex(4,3), simplex_last_1(4,3), simplex_last_2(4,3)
            !$OMP THREADPRIVATE(simplex, simplex_last_1, simplex_last_2)
            INTEGER*4, SAVE :: i, iter                    ! simplex�еĵ���                 
            LOGICAL*1, SAVE :: isOver                    ! ��ʾ�Ƿ�����ײ���߼�����
            !$OMP THREADPRIVATE(i, iter, isOver)
            REAL*8, SAVE :: dir(3), O(3), V1V2_vtr(3), V2V3_vtr(3), V3V4_vtr(3), VO_vtr(3)
            !$OMP THREADPRIVATE(dir, O, V1V2_vtr, V2V3_vtr, V3V4_vtr, VO_vtr)
            
            !---------------------------------------------------------------
            ! ��ʼ��
            !---------------------------------------------------------------
            O = .0D0
            collision_ = .FALSE.     ! ��ʼ����ײΪδ����
            colliType_ = 0
            collision_point_      = 0.D0 ! ��ײ��
            nearest_points_       = 0.D0 ! ������
            collision_normal_     = 0.D0 ! ��ײ����
            penetration_depth_    = 0.D0 ! ��͸���

            !---------------------------------------------------------------
            ! ������ײ��⣨���ΰ������㷨��
            !---------------------------------------------------------------
            CALL RoughCollisionDetection_SphericalEnvelope(p1_, p2_, collision_)
            IF ( .NOT. collision_ ) RETURN ! ���������ཻ�������GJK���
            
            !---------------------------------------------------------------
            ! ��һ��������ʼ�����塿  
            !---------------------------------------------------------------
            iter = 0
            DO WHILE(.TRUE.)
                iter = iter + 1
                ! �����������������ôֱ�ӷ��ز���ײ
                IF ( iter > 99 ) THEN
                    collision_ = .FALSE.
                    RETURN
                END IF
                
                !---------------------------------------------------------------
                ! ��1.1 ��ʼ�������� 1  �������
                !CALL RANDOM_NUMBER(dir(:))
                dir = GET_RANDOM_UNIT_VECTOR(iter)
                
                ! �����ʼ֧��ӳ��� 1����������ӵ�simplex��
                simplex(1,:) = support_mapping(p1_, p2_, dir)
                !---------------------------------------------------------------
                ! ��1.2 ��ʼ�������� 2 �� ����
                dir = - dir
            
                ! �����ʼ֧��ӳ��� 2����������ӵ�simplex��
                simplex(2,:) = support_mapping(p1_, p2_, dir)
            
                ! �����ʼ����֧��ӳ��㶼�غϣ���ô�ı��������ֱ�����غ�Ϊֹ
                IF ( ALL( DABS( simplex(1,:) - simplex(2,:) ) < 1.D-8 ) ) THEN
                    CYCLE
                ELSE
                    EXIT
                END IF
                
            END DO
            
            !---------------------------------------------------------------
            ! ��1.3 ��ʼ�������� 3 �� �߶�ָ��ԭ�㡿
            dir = VEC_PL( O, simplex(1:2,:) )
            
            ! �����ʼ֧��ӳ��� 3����������ӵ�simplex��
            simplex(3,:) = support_mapping(p1_, p2_, dir)
            
            ! �����ʱ����ӵĵ㣬��simplex��ǰ�������غϣ�˵����ԭ�㷽���Ѿ��޷�������֧��ӳ���
            ! ��ô�ɿɷ�˹����һ��������ԭ�㣬ֱ�ӷ���δ��ײ
            IF ( ALL( DABS( simplex(3, :) - simplex(1, :) ) < 1.D-8 ) .OR. &
                 ALL( DABS( simplex(3, :) - simplex(2, :) ) < 1.D-8 ) ) THEN
                collision_ = .FALSE.
                RETURN
            END IF

            !---------------------------------------------------------------
            ! ��1.4 ��ʼ�������� 4 ��ǰ��3��ӳ��֧�ֵ���ɵ������εķ���ʸ������ָ��ԭ�����ڵ���һ�ࡿ
            ! ȷ����ĵ�λ��ʸ
            V1V2_vtr = simplex(2,:) - simplex(1,:)
            V2V3_vtr = simplex(3,:) - simplex(2,:)
            
            dir = UTZVEC( CROSS_PRODUCT_3D(V1V2_vtr, V2V3_vtr) )
            
            ! ����ָ��ԭ�������   .DOT. dir 
            ! (1) = 0 ����ԭ����3��Ƭ�ڣ�˵��ԭ�������ϣ���ô�ж���ײ������
            VO_vtr = O - simplex(3,:)
            IF ( DABS(DOT_PRODUCT(VO_vtr, dir) ) < 1.D-8 ) THEN
                IF( IS_INSIDE_PF(simplex(1:3,:), O) ) THEN
                    collision_ = .TRUE. 
                    CALL EPA_solu(version_, TOL_FF_, &
                        p1_, p2_, simplex, nearest_points_, &
                        collision_normal_, collision_point_, penetration_depth_, colliType_)
                    RETURN
                END IF
            END IF
            
            ! (2) < 0 ��dir���� ; 
            IF ( DOT_PRODUCT(VO_vtr, dir) < .0D0) dir = - dir
            
            ! �����ʼ֧��ӳ��� 4����������ӵ�simplex��
            simplex(4,:) = support_mapping(p1_, p2_, dir)
            
            ! ���������ĵ����Ǹ�֧��ӳ�����ǰ3������ͬһƽ���ϣ��򷵻ز���ײ
            IF ( DABS( DIST_PF_SIGN( simplex(4,:), simplex(1:3,:) ) ) < 1.D-8 ) THEN
                collision_ = .FALSE. 
                RETURN
            END IF
            
            ! �жϵ������Ƿ����ԭ��
            ! �������ΰ���ԭ�㣬��ô��Ϊ������ײ
            IF ( isPointInSimplex( O, simplex ) ) THEN
                collision_ = .TRUE. 
                CALL EPA_solu(version_, TOL_FF_, &
                    p1_, p2_, simplex, nearest_points_, &
                    collision_normal_, collision_point_, penetration_depth_, colliType_)
                RETURN
            END IF
            
            !---------------------------------------------------------------
            ! �����������嵥������ԭ�㷽�������  
            !---------------------------------------------------------------
            ! ��ʼ�����ι������̽����������ʼ�жϹ���û�а���ԭ�㣬
            ! ��ô���濪ʼʹ��������ԭ�㷽�����
            !---------------------------------------------------------------
            simplex_last_1 = .0D0
            simplex_last_2 = .0D0
            
            iter = 0
            DO WHILE(.TRUE.)
                
                ! �ﵽ������������δ�ҵ�����ԭ��ĵ����Σ���ǰ�汾����δ��ײ
                iter = iter + 1
                IF ( iter > 50 ) THEN
                    collision_ = .FALSE. 
                    RETURN
                END IF
                
                
                ! ��¼��ǰ����������
                simplex_last_2 = simplex_last_1
                simplex_last_1 = simplex
                simplex = update_simplex_GJK( p1_, p2_, simplex )
                
                ! ���������ĵ����Ǹ�֧��ӳ�����ǰ3������ͬһƽ���ϣ��򷵻ز���ײ
                 ! ���ȱ���3�㲻����
                IF ( NORM2(CROSS_PRODUCT_3D( simplex(2,:) - simplex(1,:), simplex(3,:) - simplex(2,:) )) < 1.D-8 ) THEN
                    collision_ = .FALSE. 
                    RETURN
                ELSE
                    IF ( DABS( DIST_PF_SIGN( simplex(4,:), simplex(1:3,:) ) ) < 1.D-8 ) THEN
                        collision_ = .FALSE. 
                        RETURN
                    END IF
                END IF
                
                ! ��������ΰ���ԭ�㣨������ͱ��ϣ�����ô�ж���ײ�������������������
                IF ( isPointInSimplex(O, simplex) ) THEN
                    collision_ = .TRUE. 
                    CALL EPA_solu(version_, TOL_FF_ , &
                        p1_, p2_, simplex, nearest_points_, &
                        collision_normal_, collision_point_, penetration_depth_, colliType_)
                    RETURN
                END IF

                ! ������ֹ������ �����β��ٷ����仯�����߲��ٲ����µĵ�����
                isOver = .FALSE.

                DO i = 1, 4, 1
                    IF ( ALL( DABS( simplex(i, :) - simplex_last_1(i, :) ) < 1.D-8) .OR. &
                         ALL( DABS( simplex(i, :) - simplex_last_2(i, :) ) < 1.D-8) ) THEN
                        isOver = .TRUE.
                    ELSE
                        isOver = .FALSE.
                        EXIT
                    END IF
                END DO

                IF (isOver) THEN
                    collision_ = .FALSE.
                    RETURN
                END IF
    
            END DO
            
            RETURN
        END SUBROUTINE GJKEPA
    
        
        SUBROUTINE EPA_solu(version_, TOL_FF_,&
            p1_, p2_, simplex_, nearest_points_, &
            collision_normal_, collision_point_, penetration_depth_, collision_info_)
            IMPLICIT NONE
            INTEGER*4, INTENT(IN) :: version_
	        REAL*8, INTENT(IN) :: TOL_FF_
            REAL*8, INTENT(IN) :: p1_(:,:), p2_(:,:)        ! ͹������
            REAL*8, INTENT(IN) :: simplex_(4,3)             ! ������ԭ��ĵ�����
            REAL*8, INTENT(OUT) :: nearest_points_(2,3)     ! ������
            REAL*8, INTENT(OUT) :: collision_normal_(3)     ! ��ײ����
            REAL*8, INTENT(OUT) :: collision_point_(3)     ! ��ײ��
            REAL*8, INTENT(OUT) :: penetration_depth_       ! ��͸���
            INTEGER*4, INTENT(OUT) :: collision_info_
            !---------------------------------------------------------------
            REAL*8, SAVE, ALLOCATABLE :: polytope(:,:,:)            ! ��չ������     (��ţ���ţ���XYZ)
            REAL*8, SAVE, ALLOCATABLE :: polytope_res(:,:,:)        ! ��չ������res
            REAL*8, SAVE :: pene_depth, nml_devi(3), collision_normal_new(3)
            LOGICAL*1, SAVE :: isExpa
            INTEGER*4, SAVE :: istat, iter
            !$OMP THREADPRIVATE(polytope, polytope_res, pene_depth, nml_devi, isExpa, istat, iter, collision_normal_new) ! 
            !---------------------------------------------------------------
            isExpa = .FALSE.
            collision_point_     = 0.D0
            nearest_points_      = 0.D0
            collision_normal_    = 0.D0
            penetration_depth_   = 0.D0
            
            pene_depth = .0D0
            nml_devi = .0D0
   
            !---------------------------------------------------------------
            ! ��ʼ��
            IF ( ALLOCATED(polytope) ) DEALLOCATE(polytope)
            ALLOCATE( polytope(4,3,3), STAT = istat )     ! 4���棬ÿ����3�����㣬ÿ������3������
            IF ( ALLOCATED(polytope_res) ) DEALLOCATE(polytope_res)
            ALLOCATE( polytope_res(6,3,3), STAT = istat ) ! 6���棬ÿ����3�����㣬ÿ������3������
            
            polytope(1,1,:) = simplex_(1,:)
            polytope(1,2,:) = simplex_(2,:)
            polytope(1,3,:) = simplex_(3,:)
            
            polytope(2,1,:) = simplex_(1,:)
            polytope(2,2,:) = simplex_(3,:)
            polytope(2,3,:) = simplex_(4,:)
            
            polytope(3,1,:) = simplex_(1,:)
            polytope(3,2,:) = simplex_(2,:)
            polytope(3,3,:) = simplex_(4,:)
            
            polytope(4,1,:) = simplex_(2,:)
            polytope(4,2,:) = simplex_(3,:)
            polytope(4,3,:) = simplex_(4,:)
            
            iter = 0
            DO WHILE(.TRUE.)
                ! �ﵽ������������δ�ҵ�����ԭ��ĵ����Σ���ǰ�汾����δ��ײ
                iter = iter + 1
                IF ( iter > 99 ) THEN
                    WRITE(UNIT=6, FMT="(A)") "EPA_solu() - The current version of the EPA algorithm does not support collisions in this case." ! 
                    PAUSE
                    RETURN
                END IF
                
                ! ��չ�����壬��˳���õ�ǰpolytope�������dist_min��nml_devi
                CALL update_expandingPolytope_EPA( p1_, p2_, polytope, isExpa, polytope_res, pene_depth, nml_devi ) 
                
                IF ( isExpa == .FALSE. ) THEN
                    ! ��ֹ��������
                    penetration_depth_ = pene_depth
                    collision_normal_ = nml_devi
                    EXIT
                END IF
            
                ! Ϊ�´ε�����׼��������չ��Ķ�����������polytope��������̬�����С�������·��䣬 polytope_res����δ����״̬
                IF(ALLOCATED(polytope)) DEALLOCATE( polytope, STAT = istat )
                !ALLOCATE( polytope( SIZE(polytope_res, 1), 3, 3), STAT = istat )
            
                polytope = polytope_res
            
                IF(ALLOCATED(polytope_res)) DEALLOCATE( polytope_res, STAT = istat )

            END DO
            
            ! �ҵ���ײ����nml_devi�ʹ�͸���pene_depth������������
            nearest_points_ = get_nearest_points(p1_, p2_, collision_normal_, penetration_depth_)
            
            ! ������ײ��
            IF ( version_ == 1 ) THEN ! �汾1 
                collision_point_ = get_collisionPoint_01(p1_, p2_, collision_normal_)
            ELSE IF ( version_ == 2 ) THEN ! �汾2
                collision_point_ = get_collisionPoint_02(p1_, p2_, collision_normal_) 
            ELSE IF ( version_ == 3 ) THEN ! �汾3
                collision_point_ = get_collisionPoint_03(p1_, p2_, collision_normal_, collision_normal_new)    
                collision_normal_ = collision_normal_new
            ELSE
                WRITE(UNIT=6, FMT="(A)") "EPA_solu() - get_collisionPoint(p1_, p2_, collision_normal_)" ! 
                PAUSE
                STOP
            END IF
            
            ! ������ײ����
            collision_info_ = get_info_collisionType(p1_, p2_, collision_normal_, TOL_FF_)
            
            RETURN
        END SUBROUTINE EPA_solu
        
        !---------------------------------------------------------------
        !
        ! (1) ��ײ��Ϣ
        !
        !---------------------------------------------------------------
        FUNCTION get_info_collisionType(p1_, p2_, collision_normal_, TOL) RESULT(res_)
            IMPLICIT NONE
            REAL*8, INTENT(IN) :: p1_(:,:), p2_(:,:)        ! ͹������
            REAL*8, INTENT(IN) :: collision_normal_(3)     ! ��ײ����
	        REAL*8, INTENT(IN) :: TOL ! ��ֵ
            INTEGER*4 :: res_ ! ��ײ���ͣ�0 - δ��ײ�� 1 - ������ײ���� �� 2 - ������ײ(��ֵ�ж�)
            
            INTEGER*4, SAVE :: i, C, D ! 
            REAL*8, SAVE :: max_dot_product, dot_product_vertex, dot_product_maxIndex
            !REAL*8, PARAMETER :: TOL = 1.D0 ! ��ֵ
            !$OMP THREADPRIVATE(i, C, D, max_dot_product, dot_product_vertex, dot_product_maxIndex) ! 
            ! ����������˵��һ����������ײ����ʼ��Ϊ1
            res_ = 1
            
            !---------------------------------------------------------------
            ! ��p1��
            ! ���ҳ�����һ��֧�ŵ㣬���ɿɷ�˹�������ֵ��
            dot_product_maxIndex = - HUGE(1.0D0)
            DO i = 1, SIZE(p1_, 1)
                dot_product_vertex = DOT_PRODUCT(collision_normal_, p1_(i,:))
                IF (dot_product_vertex > dot_product_maxIndex) THEN
                    dot_product_maxIndex = dot_product_vertex
                END IF
            END DO
                
            C = 0
            DO i = 1, SIZE(p1_, 1)
                dot_product_vertex = DOT_PRODUCT(collision_normal_, p1_(i,:))
                IF (dot_product_vertex  >  dot_product_maxIndex - TOL) THEN
                    C = C + 1
                END IF
            END DO
            
            !---------------------------------------------------------------
            ! ��p2��
            ! ���ҳ�����һ��֧�ŵ㣬���ɿɷ�˹�������ֵ��
            dot_product_maxIndex = - HUGE(1.0D0)
            DO i = 1, SIZE(p2_, 1)
                dot_product_vertex = DOT_PRODUCT( - collision_normal_, p2_(i,:))
                IF (dot_product_vertex > dot_product_maxIndex) THEN
                    dot_product_maxIndex = dot_product_vertex
                END IF
            END DO
                
            ! �ٰ����е�֧�ŵ��ҳ���
            D = 0 ! �洢p2֧�ŵ�����
            DO i = 1, SIZE(p2_, 1)
                dot_product_vertex = DOT_PRODUCT( - collision_normal_, p2_(i,:))
                IF (dot_product_vertex  >  dot_product_maxIndex - TOL) THEN
                    D = D + 1
                END IF
            END DO
            
            !---------------------------------------------------------------
            ! ������ֵ�ж��Ƿ����϶�Ϊ��Ӵ�
            !---------------------------------------------------------------
            ! ���ж���ײ�����Ƿ�������ײ��
            IF ( C >= 3 .AND. D >= 3 ) res_ = 2
            
            RETURN
        END FUNCTION get_info_collisionType
        
        
        !---------------------------------------------------------------
        !
        ! (2) ��ײ��  1 - ����  2 - ͨ��   3 - ���ת��ר��
        !
        !---------------------------------------------------------------
        !---------------------------------------------------------------
        ! Ѱ����ײ�� �汾3  (p1��� - p2ת�� ר��) 
        ! a. ��ײ��һ����ת����
        ! b. ��ײ�����Ϊ��XOY���ϵ�ͶӰ���ұ�Ϊ��λ����
        ! c. ��ײ��z�����Ϊ��p1�������һ��
        FUNCTION get_collisionPoint_03(p1_, p2_, collision_normal_, collision_normal_new_) RESULT(res_)
            IMPLICIT NONE
            REAL*8, INTENT(IN) :: p1_(:,:), p2_(:,:)        ! ͹������
            REAL*8, INTENT(IN) :: collision_normal_(3)     ! ��ײ����
            REAL*8, INTENT(OUT) :: collision_normal_new_(3)     ! ������ײ����
            REAL*8 :: res_(3)
            REAL*8 :: maxDot, vertexDot
            INTEGER*4 :: i, index_support_p2 
            !Ѱ��p2��ײ�������ǰ��֧�ŵ�
            ! ����͹�����p2�� - dir_�����ϵ���Զ��
            maxDot = - HUGE(1.D0)
            index_support_p2 = 0
            DO i = 1, SIZE(p2_, 1)
                vertexDot = DOT_PRODUCT( - collision_normal_, p2_(i,:))
                IF (vertexDot  >  maxDot - 1.D-8) THEN
                    maxDot = vertexDot
                    index_support_p2 = i
                END IF
            END DO
            res_ = p2_( index_support_p2,: ) ! a. ��ײ��һ����ת����
            res_(3) = SUM(p1_(:,3)) / REAL(SIZE(p1_,1)) ! c. ��ײ��z�����Ϊ��p1�������һ��
            
            collision_normal_new_ = collision_normal_
            collision_normal_new_(3) = 0.D0
            collision_normal_new_ = collision_normal_new_ / NORM2(collision_normal_new_) ! b. ��ײ�����Ϊ��XOY���ϵ�ͶӰ���ұ�Ϊ��λ����
            RETURN
        END FUNCTION get_collisionPoint_03
        
        
        
        ! Ѱ����ײ�� �汾2 ��ͨ�ã�
        FUNCTION get_collisionPoint_02(p1_, p2_, collision_normal_) RESULT(res_)
            IMPLICIT NONE
            REAL*8, INTENT(IN) :: p1_(:,:), p2_(:,:)        ! ͹������
            REAL*8, INTENT(IN) :: collision_normal_(3)     ! ��ײ����
            REAL*8 :: res_(3)
            REAL*8, PARAMETER :: tol = 1.D-6
            TYPE(List_Array3d), SAVE :: SPT_p1, SPT_p2 ! 
            INTEGER*4, SAVE :: n1, n2
            !$OMP THREADPRIVATE(SPT_p1, SPT_p2, n1, n2) ! 
            res_ = 0.D0
            !---------------------------------------------------------------
            CALL SPT_p1%reset()
            CALL SPT_p2%reset()
            
            CALL AddAllSupports(p1_, collision_normal_, 1.D-1, SPT_p1) 
            CALL AddAllSupports(p2_, - collision_normal_, 1.D-1, SPT_p2) 
            
            n1 = SPT_p1%getSize()
            n2 = SPT_p2%getSize()
            
            IF ( n1 == 1 .AND. n2 == 1 ) THEN
                CALL case_01(SPT_p1, SPT_p2, res_)
                
            ELSE IF ( n1 == 1 .AND. n2 >= 2 ) THEN
                CALL case_02(SPT_p1, res_)    
                
            ELSE IF ( n1 >= 2 .AND. n2 == 1 ) THEN
                CALL case_02(SPT_p2, res_)  ! ע�⴫��
                
            ELSE IF ( n1 == 2 .AND. n2 == 2 ) THEN
                CALL case_03(SPT_p1, SPT_p2, res_)
                
            ELSE IF ( n1 == 2 .AND. n2 >= 3 ) THEN
                CALL case_04(SPT_p2, SPT_p1, res_)  ! ע�⴫��
                
            ELSE IF ( n1 >= 3 .AND. n2 == 2 ) THEN
                CALL case_04(SPT_p1, SPT_p2, res_)
                
            ELSE IF ( n1 >= 3 .AND. n2 >= 3 ) THEN
                CALL case_05(SPT_p1, SPT_p2, res_)
                
            ELSE
                WRITE(UNIT=6, FMT="(A)") "ERROR - FUNCTION get_collisionPoint(p1_, p2_, collision_normal_) RESULT(res_)" ! 
                READ(UNIT=5, FMT=*)   
                STOP
            END IF
            
            CALL SPT_p1%reset()
            CALL SPT_p2%reset()

            RETURN
            CONTAINS
                SUBROUTINE AddAllSupports(p_, nml_, tol_, SPT_p)
                    IMPLICIT NONE
                    REAL*8, INTENT(IN) :: p_(:,:), nml_(3), tol_  ! ͹������
                    TYPE(List_Array3d), INTENT(OUT) :: SPT_p
                    INTEGER*4, SAVE :: i
                    REAL*8, SAVE :: dot_maxLoc, dot_now
                    !$OMP THREADPRIVATE(i, dot_maxLoc, dot_now) ! 
                    ! ���ҳ�����һ��֧�ŵ㣬���ɿɷ�˹�������ֵ��
                    dot_maxLoc = - HUGE(1.0D0)
                    DO i = 1, SIZE(p_, 1)
                        dot_now = DOT_PRODUCT(nml_, p_(i,:))
                        IF (dot_now > dot_maxLoc) dot_maxLoc = dot_now
                    END DO
                    ! �ٰ����е�֧�ŵ��ҳ���
                    CALL SPT_p%reset()
                    DO i = 1, SIZE(p_, 1)
                        dot_now = DOT_PRODUCT(nml_, p_(i,:))
                        IF (dot_now  >  dot_maxLoc - tol_) CALL SPT_p%append( p_(i,:) )
                    END DO
                    RETURN
                END SUBROUTINE AddAllSupports
        
                SUBROUTINE case_01(SPT_p1_, SPT_p2_, colliPoin_) ! ��ײ��Ϊ�������ߵ��е�
                    IMPLICIT NONE
                    TYPE(List_Array3d), INTENT(IN) :: SPT_p1_, SPT_p2_ ! 
                    REAL*8, INTENT(OUT) :: colliPoin_(3)    
                    TYPE(ListNode_Array3d), SAVE :: node(2) ! 
                    !$OMP THREADPRIVATE(node) ! 
                    node(1) = SPT_p1_%atIndex(1)
                    node(2) = SPT_p2_%atIndex(1)
                    colliPoin_ = ( node(1)%value(:) + node(2)%value(:) ) / 2.0D0
                    RETURN
                END SUBROUTINE case_01
                
                SUBROUTINE case_02(SPT_p_, colliPoin_) ! ��ײ���� p1 ��
                    IMPLICIT NONE
                    TYPE(List_Array3d), INTENT(IN) :: SPT_p_ ! 
                    REAL*8, INTENT(OUT) :: colliPoin_(3)    
                    TYPE(ListNode_Array3d), SAVE :: node ! 
                    !$OMP THREADPRIVATE(node) ! 
                    node = SPT_p_%atIndex(1)
                    colliPoin_ = node%value(:)
                    RETURN
                END SUBROUTINE case_02
                
                SUBROUTINE case_03(SPT_p1_, SPT_p2_, colliPoin_) ! ��������
                    IMPLICIT NONE
                    TYPE(List_Array3d), INTENT(IN) :: SPT_p1_, SPT_p2_ ! 
                    REAL*8, INTENT(OUT) :: colliPoin_(3)    
                    TYPE(ListNode_Array3d), SAVE :: node1(2), node2(2) !    
                    REAL*8, SAVE :: foot(2,3), ep1(2,3), ep2(2,3)
                    !$OMP THREADPRIVATE(node1, node2, foot, ep1, ep2) ! 
                    node1(1) = SPT_p1_%atIndex(1)
                    node1(2) = SPT_p1_%atIndex(2)
                    node2(1) = SPT_p2_%atIndex(1)
                    node2(2) = SPT_p2_%atIndex(2)
                    ep1(1,:) = node1(1)%value(:)
                    ep1(2,:) = node1(2)%value(:)
                    ep2(1,:) = node2(1)%value(:)
                    ep2(2,:) = node2(2)%value(:)
                    foot = FOOT_LL(ep1, ep2)
                    colliPoin_ = ( foot(1,:) + foot(2,:) ) / 2.0D0
                    RETURN
                END SUBROUTINE case_03
                
                
                SUBROUTINE case_04(SPT_p1_, SPT_p2_, colliPoin_)
                    IMPLICIT NONE
                    TYPE(List_Array3d), INTENT(IN) :: SPT_p1_, SPT_p2_ ! >=3   =2  
                    REAL*8, INTENT(OUT) :: colliPoin_(3)
                    TYPE(ListNode_Array3d), SAVE, ALLOCATABLE :: node1(:) ! 
                    TYPE(ListNode_Array3d), SAVE :: node2(2) ! 
                    INTEGER*4, SAVE :: i, sz
                    REAL*8, SAVE, ALLOCATABLE :: sprt1(:,:)
                    REAL*8, SAVE :: XYZcent(3), sprt2(2,3)
                    !$OMP THREADPRIVATE(node1, node2, i, sz, sprt1, XYZcent, sprt2) ! 
                    sz = SPT_p1_%getSize()
                    IF(ALLOCATED(node1)) DEALLOCATE(node1)
                    ALLOCATE( node1( sz ) )
                    IF(ALLOCATED(sprt1)) DEALLOCATE(sprt1)
                    ALLOCATE( sprt1( sz, 3 ) )
                    DO i = 1, sz
                        node1(i) = SPT_p1_%atIndex(i)
                        sprt1(i,:) = node1(i)%value(:)
                    END DO
                    node2(1) = SPT_p2_%atIndex(1)
                    node2(2) = SPT_p2_%atIndex(2)
                    sprt2(1,:) = node2(1)%value(:)
                    sprt2(2,:) = node2(2)%value(:)
                    
                    SELECT CASE( branch_case_04(sprt1, sprt2)  )
                        CASE (1)
                            CALL case_04_1(sprt1, sprt2, colliPoin_)
                        CASE (2)
                            CALL case_04_2(sprt1, sprt2, colliPoin_)
                        CASE (3)
                            CALL case_04_3(sprt1, sprt2, colliPoin_) ! �����£�Ŀǰ��case_04_2ִ��ͬ���Ĳ���
                        CASE DEFAULT
                        
                    END SELECT

                    DEALLOCATE(node1)
                    RETURN
                END SUBROUTINE case_04
                    
                    FUNCTION branch_case_04(sprt1, sprt2) RESULT(res_)
                        IMPLICIT NONE
                        REAL*8, INTENT(IN) :: sprt1(:,:), sprt2(2,3) ! >=3   =2 
                        INTEGER*4 :: res_
                        REAL*8, SAVE, ALLOCATABLE :: sprt1_sorted(:,:)
                        INTEGER*4, SAVE :: i, C ! 
                        !$OMP THREADPRIVATE(sprt1_sorted, i, C ) ! 
                        IF(ALLOCATED(sprt1_sorted)) DEALLOCATE(sprt1_sorted)
                        sprt1_sorted = SORT_CLOCK(sprt1)
                        C = 0
                        DO i = 1, 2 
                            IF( IS_INSIDE_PF(sprt1_sorted, sprt2(i,:)) ) C = C + 1
                        END DO
                        SELECT CASE( C )
                            CASE (0)
                                res_ = 1
                            CASE (1)
                                res_ = 3
                            CASE (2)
                                res_ = 2
                            CASE DEFAULT
                                WRITE(UNIT=6, FMT="(A)") "ERROR - FUNCTION branch_case_04(sprt1, sprt2) RESULT(res_)" ! 
                                READ(UNIT=5, FMT=*)   
                                STOP
                        END SELECT
                        RETURN
                    END FUNCTION branch_case_04
                
                    SUBROUTINE case_04_1(sprt1, sprt2, colliPoin_)
                        IMPLICIT NONE
                        REAL*8, INTENT(IN) :: sprt1(:,:), sprt2(2,3) ! >=3   =2 
                        REAL*8, INTENT(OUT) :: colliPoin_(3)    
                        INTEGER*4, SAVE :: i, sz
                        REAL*8, SAVE :: XYZcent(3)
                        !$OMP THREADPRIVATE(i, sz, XYZcent) ! 
                        sz = SIZE(sprt1, 1)
                        FORALL(i = 1:3) XYZcent(i) = SUM( sprt1(:,i) ) / REAL(sz, KIND(1.D0)) 
                        colliPoin_ = FOOT_PL(XYZcent, sprt2)
                        RETURN
                    END SUBROUTINE case_04_1
                        
                    SUBROUTINE case_04_2(sprt1, sprt2, colliPoin_)
                        IMPLICIT NONE
                        REAL*8, INTENT(IN) :: sprt1(:,:), sprt2(2,3) ! >=3   =2 
                        REAL*8, INTENT(OUT) :: colliPoin_(3)    
                        colliPoin_ = (sprt2(1,:) + sprt2(2,:)) * .5D0
                        RETURN
                    END SUBROUTINE case_04_2
                        
                    SUBROUTINE case_04_3(sprt1, sprt2, colliPoin_)
                        IMPLICIT NONE
                        REAL*8, INTENT(IN) :: sprt1(:,:), sprt2(2,3) ! >=3   =2 
                        REAL*8, INTENT(OUT) :: colliPoin_(3)    
                        colliPoin_ = (sprt2(1,:) + sprt2(2,:)) * .5D0
                        RETURN
                    END SUBROUTINE case_04_3
                
                
                SUBROUTINE case_05(SPT_p1_, SPT_p2_, colliPoin_) ! �����Ż�����ǰ���߼�����ײ����p1�ϣ�Ϊ����sprt����ƽ��
                    IMPLICIT NONE
                    TYPE(List_Array3d), INTENT(IN) :: SPT_p1_, SPT_p2_ ! 
                    REAL*8, INTENT(OUT) :: colliPoin_(3)    
                    TYPE(ListNode_Array3d), SAVE, ALLOCATABLE :: node1(:) ! 
                    TYPE(ListNode_Array3d), SAVE :: node2(2) ! 
                    INTEGER*4, SAVE :: i, sz
                    REAL*8, SAVE, ALLOCATABLE :: sprt1(:,:)
                    REAL*8, SAVE :: sprt2(2,3)
                    !$OMP THREADPRIVATE(node1, node2, i, sz, sprt1, sprt2) ! 
                    sz = SPT_p1_%getSize()
                    IF(ALLOCATED(node1)) DEALLOCATE(node1)
                    ALLOCATE( node1( sz ) )
                    IF(ALLOCATED(sprt1)) DEALLOCATE(sprt1)
                    ALLOCATE( sprt1( sz, 3 ) )
                    DO i = 1, sz
                        node1(i) = SPT_p1_%atIndex(i)
                        sprt1(i,:) = node1(i)%value(:)
                    END DO
                    FORALL(i = 1:3) colliPoin_(i) = SUM( sprt1(:,i) ) / REAL(sz, KIND(1.D0)) 
                    DEALLOCATE(node1)
                    RETURN
                END SUBROUTINE case_05
                
        END FUNCTION get_collisionPoint_02

        !---------------------------------------------------------------
        ! Ѱ����ײ�� �汾1 �����԰棩
        FUNCTION get_collisionPoint_01(p1_, p2_, collision_normal_) RESULT(res_)
            IMPLICIT NONE
            REAL*8, INTENT(IN) :: p1_(:,:), p2_(:,:)        ! ͹������
            REAL*8, INTENT(IN) :: collision_normal_(3)     ! ��ײ����
            REAL*8 :: res_(3)

            INTEGER*4, SAVE :: index_support_p1(2), index_support_p2(2) ! ������������������֧�ŵ�����
            INTEGER*4, SAVE :: i, istat, C, D, max_index 
            REAL*8, SAVE :: maxDot, vertexDot, maxDotIndex
            REAL*8, SAVE :: foot_points(2,3), lineDefiEP_1(2,3), lineDefiEP_2(2,3)
            REAL*8, SAVE, ALLOCATABLE :: supports_for_aver(:,:)
            !$OMP THREADPRIVATE(index_support_p1, index_support_p2, i, istat, C, D, max_index ) ! 
            !$OMP THREADPRIVATE(maxDot, vertexDot, maxDotIndex ) ! 
            !$OMP THREADPRIVATE(foot_points, lineDefiEP_1, lineDefiEP_2 ) ! 
            !$OMP THREADPRIVATE(supports_for_aver ) ! 
            
            IF(ALLOCATED(supports_for_aver)) DEALLOCATE(supports_for_aver)
            
            res_ = 0.D0
            
            !Ѱ��p1��ײ�������ǰ��2��֧�ŵ�
            ! ����͹�����p1��dir_�����ϵ���Զ��
            maxDot = - HUGE(1.D0)
            index_support_p1 = 0
            DO i = 1, SIZE(p1_, 1)
                vertexDot = DOT_PRODUCT(collision_normal_, p1_(i,:))
                IF (vertexDot  >  maxDot - 1.D-8) THEN
                    maxDot = vertexDot
                    index_support_p1(2) =  index_support_p1(1)
                    index_support_p1(1) = i
                END IF
            END DO
            IF( index_support_p1(2) == 0 ) index_support_p1(2) = index_support_p1(1) ! �����һ�ξ��ҵ������ֵ����ô��Ҫ������=0
            
            
            !Ѱ��p2��ײ�������ǰ��2��֧�ŵ�
            ! ����͹�����p2��-dir_�����ϵ���Զ��
            maxDot = - HUGE(1.D0)
            index_support_p2 = 0
            DO i = 1, SIZE(p2_, 1)
                vertexDot = DOT_PRODUCT( - collision_normal_, p2_(i,:))
                IF (vertexDot  >  maxDot - 1.D-8) THEN
                    maxDot = vertexDot
                    index_support_p2(2) =  index_support_p2(1)
                    index_support_p2(1) = i
                END IF
            END DO
            IF( index_support_p2(2) == 0 ) index_support_p2(2) = index_support_p2(1) ! �����һ�ξ��ҵ������ֵ����ô��Ҫ������=0
            
            
            !Ѱ��p2��ײ����������֧�ŵ�
            
            
            ! ��case 1 ������p��ֻ�ѵ���һ��֧�ŵ㣬��ȡ�������е�Ϊ��ײ��
            IF ( index_support_p1(1) == index_support_p1(2) .AND. index_support_p2(1) == index_support_p2(2) ) THEN
                res_ = ( p1_( index_support_p1(1), :) + p2_( index_support_p2(1), :) ) / 2.D0
            END IF
            
            ! ��case 2 ��һ���ѵ��˶��֧�ŵ㣬����һ��convex����ײ������ֻ��������һ��֧�ŵ㣬��õ���Ϊ��ײ��
            IF ( index_support_p1(1) /= index_support_p1(2) .AND. index_support_p2(1) == index_support_p2(2) ) THEN
                res_ = p2_( index_support_p2(1), :)
            ELSE IF (  index_support_p1(1) == index_support_p1(2) .AND. index_support_p2(1) /= index_support_p2(2) ) THEN
                res_ = p1_( index_support_p1(1), :)
            END IF
            
            ! ��case 3 ������convex����ײ�����϶��ѵ��˶��֧�ŵ㣬���ڸ���������϶���ײ����p1�ϣ����ݶ�����ײ��Ϊ��ײ����������֧�ŵ������ƽ����
            IF ( index_support_p1(1) /= index_support_p1(2) .AND. index_support_p2(1) /= index_support_p2(2) ) THEN
                !---------------------------------------------------------------
                !Ѱ��p1��ײ�������ǰ��2��֧�ŵ�
                ! ����͹�����p1��dir_�����ϵ���Զ��
                maxDot = - HUGE(1.D0)
                index_support_p1 = 0
                C = 0 ! �洢p1֧�ŵ�����
                ALLOCATE( supports_for_aver( SIZE(p1_, 1), 3 ), STAT = istat )
                supports_for_aver = 0.D0
                
                !---------------------------------------------------------------
                ! ��p1��
                ! ���ҳ�����һ��֧�ŵ㣬���ɿɷ�˹�������ֵ��
                maxDotIndex = - HUGE(1.0D0)
                max_index = 1
                DO i = 1, SIZE(p1_, 1)
                    vertexDot = DOT_PRODUCT(collision_normal_, p1_(i,:))
                    IF (vertexDot > maxDotIndex) THEN
                        maxDotIndex = vertexDot
                        max_index = i
                    END IF
                END DO
                
                ! �ٰ����е�֧�ŵ��ҳ���
                DO i = 1, SIZE(p1_, 1)
                    vertexDot = DOT_PRODUCT(collision_normal_, p1_(i,:))
                    IF (vertexDot  >  maxDotIndex - 1.D-1) THEN
                        C = C + 1
                        supports_for_aver(C,:) = p1_(i,:)
                    END IF
                END DO
        
                !���ҳ����ĵ�ȡ��ƽ�������Ż����˴�Ӧ������ײ����ཻ������ġ�
                FORALL(i = 1:3) res_(i) = SUM( supports_for_aver(1:C,i) ) / REAL(C, KIND(1.D0)) 
                
                ! �ͷſռ�
                DEALLOCATE(supports_for_aver)
                
            END IF
            RETURN
        END FUNCTION get_collisionPoint_01
        
        !---------------------------------------------------------------
        !
        ! (2) ������
        !
        !---------------------------------------------------------------
        FUNCTION get_nearest_points(p1_, p2_, collision_normal_, penetration_depth_) RESULT(res_)
        IMPLICIT NONE
            REAL*8 :: res_(2,3)
            REAL*8, INTENT(IN) :: p1_(:,:), p2_(:,:)        ! ͹������
            REAL*8, INTENT(IN) :: collision_normal_(3)     ! ��ײ����
            REAL*8, INTENT(IN) :: penetration_depth_       ! ��͸���
            !---------------------------------------------------------------
            REAL*8, SAVE :: dir_(3), support(2,3)
            INTEGER*4, SAVE :: i, max_index1, max_index2               ! ��ʱ������ѭ������������������
            REAL*8, SAVE :: max_dot_product1, max_dot_product2, dot_product_temp  ! ��ʱ���������ֵ
            !$OMP THREADPRIVATE(dir_, support, i, max_index1, max_index2, max_dot_product1, max_dot_product2, dot_product_temp) ! 
            dir_ = collision_normal_
            
            ! ����͹�����p1��dir_�����ϵ���Զ��
            max_dot_product1 = - HUGE(1.0D0)
            max_index1 = 1
            DO i = 1, SIZE(p1_, 1)
                dot_product_temp = DOT_PRODUCT(dir_, p1_(i,:))
                IF (dot_product_temp > max_dot_product1) THEN
                    max_dot_product1 = dot_product_temp
                    max_index1 = i
                END IF
            END DO
            
            ! ����͹�����p2��-dir_�����ϵ���Զ��
            max_dot_product2 = - HUGE(1.0D0)
            max_index2 = 1
            DO i = 1, SIZE(p2_, 1)
                dot_product_temp = DOT_PRODUCT(-dir_, p2_(i,:))
                IF (dot_product_temp > max_dot_product2) THEN
                    max_dot_product2 = dot_product_temp
                    max_index2 = i
                END IF
            END DO
            
            support(1,:) = p1_(max_index1,:)
            support(2,:) = p2_(max_index2,:)
            
            res_(1,:) = support(1,:) !- .5D0 * penetration_depth_ * collision_normal_
            res_(2,:) = support(2,:) !+ .5D0 * penetration_depth_ * collision_normal_
            
            RETURN
        END FUNCTION get_nearest_points
        
    
        !---------------------------------------------------------------
        !
        !  ����ԭ�㷽����չ������
        !  
        !---------------------------------------------------------------
        SUBROUTINE update_expandingPolytope_EPA( p1_, p2_, polytope_1_, isExpa_, polytope_2_, penetration_depth_, normDeviOrig_ ) 
            IMPLICIT NONE
            REAL*8, INTENT(IN) :: p1_(:,:), p2_(:,:)                              ! �������������͹����εĶ��㼯
            REAL*8, INTENT(IN) :: polytope_1_(:,:,:)                              ! ��չǰ�Ķ�����
            LOGICAL*1, INTENT(OUT) :: isExpa_                                     ! ��������Ķ�����ĳ�湲�棬�򷵻�false
            REAL*8, INTENT(OUT), ALLOCATABLE :: polytope_2_(:,:,:)  ! ��չ��Ķ�����
            REAL*8, INTENT(OUT) :: penetration_depth_, normDeviOrig_(3)                  ! ���ԭ�㵽��ľ��룬����ԭ��ķ���ʸ��
            !---------------------------------------------------------------
            REAL*8, PARAMETER :: O(3) = [0.D0, 0.D0, 0.D0]
            REAL*8, SAVE :: SPMP(3), dir(3)
            INTEGER*4, SAVE :: i, j, istat, info, n ! 
            REAL*8, SAVE, ALLOCATABLE :: dist1(:), dist2(:)
            INTEGER*4, SAVE :: min_loc_1(1), min_loc_2(1)
            REAL*8, SAVE :: min_val_1, min_val_2, dot, M(3)
            REAL*8, SAVE :: face_dele(3,3), temp
            REAL*8, SAVE, ALLOCATABLE :: scatPoints(:,:), scatPoints_temp(:,:)
            
            !$OMP THREADPRIVATE(SPMP, dir, i, j, istat, info, n, dist1, dist2, min_loc_1, min_loc_2 ) ! 
            !$OMP THREADPRIVATE(min_val_1, min_val_2, dot, M, face_dele, temp, scatPoints, scatPoints_temp ) ! 
            !---------------------------------------------------------------
            ! ��ʼ��
            isExpa_ = .FALSE.
            penetration_depth_ = 0.D0
            normDeviOrig_ = 0.D0

            IF( ALLOCATED(dist1) ) DEALLOCATE(dist1)
            ALLOCATE( dist1( SIZE( polytope_1_, 1) ), STAT = istat )
            !---------------------------------------------------------------
            ! �õ���ǰ�������չ��������浽ԭ��Ĵ�ֱ����

            DO i = 1, SIZE( polytope_1_, 1), 1
                dist1(i) = DABS( DIST_PF_SIGN( O, polytope_1_(i, :, :) ) )
            END DO
            
            ! �ҵ���̾����Ӧ���Ǹ���
            min_loc_1 = MINLOC( dist1 )
            min_val_1 = MINVAL( dist1 )
            
            ! �����Ӧ����汳��ԭ��ĵ�λ����ʸ��
            dir = UNINML( polytope_1_( min_loc_1(1), :, :) )
            dot = DOT_PRODUCT( polytope_1_( min_loc_1(1), 1, :) - O, dir)
              ! ע�⣬���ԭ����������ı����ϣ���ô�������ĵ���Ϊ�жϷ���Ļ�׼�㣬����dot=0ʱ����������ķ����෴���޷�����⵽
            IF ( DABS( dot ) < 1.D-12 ) THEN
                FORALL(i = 1:3)  M(i) = SUM( polytope_1_(:,:,i) ) / ( SIZE(polytope_1_, 1) * SIZE(polytope_1_, 2) )
                dot = DOT_PRODUCT( polytope_1_( min_loc_1(1), 1, :) - M, dir)
            END IF
                ! ���dotС��0����ôdir����
            IF( dot <= - 1.D-12 ) dir = - dir
                
            
            ! ����ʸ��������Ѱ֧��ӳ���
            SPMP = support_mapping(p1_, p2_, dir)
            
            
            !---------------------------------------------------------------
            ! ������ -> ɢ��
            IF( ALLOCATED(scatPoints) ) DEALLOCATE(scatPoints)
            CALL getHullMeshesVertex(polytope_1_, scatPoints, info)
            
            ! ֧��ӳ�����ӽ�ɢ���飬����ɢ����
            IF( ALLOCATED(scatPoints_temp) ) DEALLOCATE(scatPoints_temp)
            scatPoints_temp = scatPoints
            IF( ALLOCATED(scatPoints) ) DEALLOCATE(scatPoints)
            
            ALLOCATE( scatPoints( 1 + SIZE(scatPoints_temp, 1), 3 ), STAT = istat )
            scatPoints( 1 : SIZE(scatPoints_temp, 1), :) = scatPoints_temp(:,:)
            
            !CALL MOVE_ALLOC(scatPoints_temp, scatPoints)
            
            scatPoints( 1 + SIZE(scatPoints_temp, 1), :) = SPMP(:)
            
            !! ���dist == 0 ��ô˵��ԭ�������ϣ��޷��ж���ײ����ָ����һ�࣬��ʱ��������������һ��֧�ŵ����ɢ����
            IF ( DABS(min_val_1) < 1.D-12 ) THEN
                SPMP = support_mapping(p1_, p2_, - dir)
                ! ֧��ӳ�����ӽ�ɢ���飬����ɢ����
                IF( ALLOCATED(scatPoints_temp) ) DEALLOCATE(scatPoints_temp)
                scatPoints_temp = scatPoints
                IF( ALLOCATED(scatPoints) ) DEALLOCATE(scatPoints)
                ALLOCATE( scatPoints( 1 + SIZE(scatPoints_temp, 1), 3 ), STAT = istat )
                scatPoints( 1 : SIZE(scatPoints_temp, 1), :) = scatPoints_temp(:,:)
                scatPoints( 1 + SIZE(scatPoints_temp, 1), :) = SPMP(:)
            END IF
            
            

            ! ����͹������
            IF( ALLOCATED(polytope_2_) ) DEALLOCATE(polytope_2_)
            CALL QuickHull(scatPoints, polytope_2_, info)

            
            
            !---------------------------------------------------------------
            ! �õ���չ��Ķ�������浽ԭ��Ĵ�ֱ����
            IF( ALLOCATED(dist2) ) DEALLOCATE(dist2)
            ALLOCATE( dist2( SIZE( polytope_2_, 1) ), STAT = istat )
            DO i = 1, SIZE( polytope_2_, 1), 1
                dist2(i) = DABS( DIST_PF_SIGN( O, polytope_2_(i, :, :) ) )
            END DO
            
            ! �ҵ���̾����Ӧ���Ǹ���
            min_loc_2 = MINLOC( dist2 )
            min_val_2 = MINVAL( dist2 )
            
            ! �����Ӧ����汳��ԭ��ĵ�λ����ʸ��
            dir = UNINML( polytope_2_( min_loc_2(1), :, :) )
            dot = DOT_PRODUCT( polytope_2_( min_loc_2(1), 1, :) - O, dir)
            IF( dot < 0.D0 ) dir = - dir
            
            ! �ж��Ƿ�����չ.(ע�⣺�������Ҫ����ʽ����������ײ�����ɿɷ�˹�������ĵ���ֹͣ��������ĵ�)
            IF ( (SIZE(dist1, 1) == SIZE(dist2, 1)) ) THEN
                n = SIZE(dist1, 1)
                ! ��dist1��dist2ð������
                DO i = 1, n - 1
                    DO j = 1, n - i
                        IF (dist1(j) > dist1(j + 1)) THEN
                            ! ����Ԫ��
                            temp = dist1(j)
                            dist1(j) = dist1(j + 1)
                            dist1(j + 1) = temp
                        END IF
                        
                        IF (dist2(j) > dist2(j + 1)) THEN
                            ! ����Ԫ��
                            temp = dist2(j)
                            dist2(j) = dist2(j + 1)
                            dist2(j + 1) = temp
                        END IF
                    END DO
                END DO
                
                ! �жϴ洢�������������������Ƿ�δ�ڱ仯
                IF ( ALL( DABS(dist1 - dist2) < 1.D-8) ) THEN
                    isExpa_ = .FALSE. ! �����ⲿ���Ѿ���������չ�ˣ���ǰ�Ľ�����Ⱦ��Ƿ������ߵ���С����
                    ! ��������
                    penetration_depth_ = min_val_2
                    normDeviOrig_ = dir
                ELSE
                    isExpa_ = .TRUE.  ! �����ⲿ�����ܼ�����չ�����Լ�������
                    penetration_depth_ = 0.D0
                    normDeviOrig_ = 0.D0
                END IF
                
            ELSE IF ( (SIZE(dist1, 1) > SIZE(dist2, 1)) ) THEN ! ˵��QuickHull�̵���һ���ر𿿽����ϵ㣬ֱ�ӷ���
                isExpa_ = .FALSE. ! �����ⲿ���Ѿ���������չ�ˣ���ǰ�Ľ�����Ⱦ��Ƿ������ߵ���С����
                ! ��������
                penetration_depth_ = min_val_2
                normDeviOrig_ = dir

            ELSE
                isExpa_ = .TRUE.  ! �����ⲿ�����ܼ�����չ�����Լ�������
                penetration_depth_ = 0.D0
                normDeviOrig_ = 0.D0
            END IF

            ! �ͷ��ڴ�
            IF(ALLOCATED(dist1)) DEALLOCATE( dist1 )
            IF(ALLOCATED(dist2)) DEALLOCATE( dist2 )
            
            RETURN
        END SUBROUTINE update_expandingPolytope_EPA
        
        
        !---------------------------------------------------------------
        !
        ! ����֧��ӳ��㣨�ɿɷ�˹���
        !
        !---------------------------------------------------------------
        FUNCTION support_mapping(p1_, p2_, dir_) RESULT(res_)
            IMPLICIT NONE
            REAL(8), INTENT(IN) :: p1_(:,:), p2_(:,:), dir_(3)  ! �������������͹����εĶ��㼯����������
            REAL(8) :: res_(3)                                 ! ������������õ���֧��ӳ���
            INTEGER, SAVE :: i, maxIndex1, maxIndex2               ! ��ʱ������ѭ������������������
            REAL(8), SAVE :: maxDot1, maxDot2, tempDot  ! ��ʱ���������ֵ
            !$OMP THREADPRIVATE(i, maxIndex1, maxIndex2, maxDot1, maxDot2, tempDot ) ! 
            ! ����͹�����p1��dir_�����ϵ���Զ��
            maxDot1 = - HUGE(1.0D0)
            maxIndex1 = 1
            DO i = 1, SIZE(p1_, 1)
                tempDot = DOT_PRODUCT(dir_, p1_(i,:))
                IF (tempDot > maxDot1) THEN
                    maxDot1 = tempDot
                    maxIndex1 = i
                END IF
            END DO
            
            ! ����͹�����p2��-dir_�����ϵ���Զ��
            maxDot2 = - HUGE(1.0D0)
            maxIndex2 = 1
            DO i = 1, SIZE(p2_, 1)
                tempDot = DOT_PRODUCT(-dir_, p2_(i,:))
                IF (tempDot > maxDot2) THEN
                    maxDot2 = tempDot
                    maxIndex2 = i
                END IF
            END DO
            
            ! ����֧��ӳ���
            res_ = p1_(maxIndex1,:) - p2_(maxIndex2,:)

        END FUNCTION support_mapping

        
        !---------------------------------------------------------------
        !
        !  ��ԭ�㷽����µ����Σ������壩�����µ��Ǹ�����simplex_(4,:)������3�������ĵ���simplex_(1:3,:)
        !
        !---------------------------------------------------------------
        FUNCTION update_simplex_GJK( p1_, p2_, simplex_ ) RESULT(res_)
            IMPLICIT NONE
            REAL*8, INTENT(IN) :: p1_(:,:), p2_(:,:)  ! �������������͹����εĶ��㼯
            REAL*8, INTENT(IN) :: simplex_(4,3)  ! 
            REAL*8 :: res_(4,3)     ! ������������º�ĵ�����
            REAL*8 :: dir(3)                     ! ����
            REAL*8 :: M(3)                     ! �����ε�����ƽ��
            REAL*8 :: MO(3)                     ! �����ε�����ƽ�� -> ԭ��
            REAL*8 :: O(3) !ԭ��
            REAL*8 :: AB(3), BC(3), nml(4,3)  ! ��ʱ�����������ͷ�����
            REAL*8 :: dist_with_sign_total(4)
            REAL*8 :: SM(3)
            INTEGER*4 :: i ! 
            INTEGER*4 :: max_location(1)
            !---------------------------------------------------------------
            ! �����ε�����ƽ���㣬�����б���
            FORALL(i = 1:3) M(i) = SUM( simplex_(:,i) ) / 4.D0            
            MO = - M
            O = [0.D0, 0.D0, 0.D0]
            
            ! ����������4���泯��ĵ�λ��ʸ������->ԭ��ʸ���ĵ��
            !�泯��ķ�ʸ
            ! face 1
            AB = simplex_(1,:) - simplex_(3,:)
            BC = simplex_(3,:) - simplex_(4,:)
            nml(1,:) = UTZVEC( CROSS_PRODUCT_3D( AB, BC ) )
            IF ( DOT_PRODUCT( nml(1,:), simplex_(1,:) - M  ) < .0D0 ) nml(1,:) = - nml(1,:)
            dist_with_sign_total(1) = DOT_PRODUCT( - nml(1,:), simplex_(1,:) - O  )
            
            ! face 2
            AB = simplex_(1,:) - simplex_(2,:)
            BC = simplex_(2,:) - simplex_(4,:)
            nml(2,:) = UTZVEC( CROSS_PRODUCT_3D( AB, BC ) )
            IF ( DOT_PRODUCT( nml(2,:), simplex_(1,:) - M ) < .0D0 ) nml(2,:) = - nml(2,:)
            dist_with_sign_total(2) = DOT_PRODUCT( - nml(2,:), simplex_(1,:) - O  )
            
            ! face 3
            AB = simplex_(1,:) - simplex_(2,:)
            BC = simplex_(2,:) - simplex_(3,:)
            nml(3,:) = UTZVEC( CROSS_PRODUCT_3D( AB, BC ) )
            IF ( DOT_PRODUCT( nml(3,:), simplex_(1,:) - M ) < .0D0 ) nml(3,:) = - nml(3,:)
            dist_with_sign_total(3) = DOT_PRODUCT( - nml(3,:), simplex_(1,:) - O  )
            
            ! face 4
            AB = simplex_(2,:) - simplex_(3,:)
            BC = simplex_(3,:) - simplex_(4,:)
            nml(4,:) = UTZVEC( CROSS_PRODUCT_3D( AB, BC ) )
            IF ( DOT_PRODUCT( nml(4,:), simplex_(2,:) - M ) < .0D0 ) nml(4,:) = - nml(4,:)
            dist_with_sign_total(4) = DOT_PRODUCT( - nml(4,:), simplex_(2,:) - O  )
            
            ! �ҳ�������ԭ����Ǹ���
            max_location = MAXLOC( dist_with_sign_total )
            
            ! �жϳ���һ����������
            dir = nml( max_location(1), :)
            
            ! ��һ��֧��ӳ���
            SM = support_mapping(p1_, p2_, dir)
            
            ! ���ص�����
            SELECT CASE( max_location(1) )
              CASE (1)
                res_(1,:) = simplex_(1,:)
                res_(2,:) = simplex_(3,:)
                res_(3,:) = simplex_(4,:)
                res_(4,:) = SM
              CASE (2)
                res_(1,:) = simplex_(1,:)
                res_(2,:) = simplex_(2,:)
                res_(3,:) = simplex_(4,:)
                res_(4,:) = SM
              CASE (3)
                res_(1,:) = simplex_(1,:)
                res_(2,:) = simplex_(2,:)
                res_(3,:) = simplex_(3,:)
                res_(4,:) = SM
              CASE (4)
                res_(1,:) = simplex_(2,:)
                res_(2,:) = simplex_(3,:)
                res_(3,:) = simplex_(4,:)
                res_(4,:) = SM
              CASE DEFAULT
                
              END SELECT
            
              RETURN
            
        END FUNCTION update_simplex_GJK
        
        
        !---------------------------------------------------------------
        !
        ! �����������ײ���
        !
        !---------------------------------------------------------------
        SUBROUTINE RoughCollisionDetection_SphericalEnvelope(p1_, p2_, isColl_)
            IMPLICIT NONE
            REAL*8, INTENT(IN) :: p1_(:,:), p2_(:,:)  ! �������������͹����εĶ��㼯
            LOGICAL*1, INTENT(OUT) :: isColl_ ! �Ƿ���ײ
            REAL*8, SAVE :: mp1(3), mp2(3), r1, r2
            REAL*8 :: dist1(SIZE(p1_,1)), dist2(SIZE(p2_,1))
            INTEGER*4, SAVE :: i ! 
            REAL*8, PARAMETER :: TOL = 1.D0
            !$OMP THREADPRIVATE(mp1, mp2, r1, r2, i) ! 
            ! Ѱ������ƽ����
            FORALL (i = 1:3) mp1(i) = SUM( p1_(:,i) ) / SIZE(p1_, 1)
            FORALL (i = 1:3) mp2(i) = SUM( p2_(:,i) ) / SIZE(p2_, 1)
            
            ! Ѱ������뾶
            FORALL(i = 1:SIZE(p1_,1)) dist1(i) = NORM2(p1_(i,:) - mp1(:))
            r1 = MAXVAL( dist1 )
            FORALL(i = 1:SIZE(p2_,1)) dist2(i) = NORM2(p2_(i,:) - mp2(:))
            r2 = MAXVAL( dist2 )
            
            ! �ж������Ƿ��ཻ
            isColl_ = MERGE(.TRUE., .FALSE., NORM2(mp1 - mp2) <= r1 + r2 + TOL)
            
            RETURN
        END SUBROUTINE RoughCollisionDetection_SphericalEnvelope

        
        
        
        !_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
        !_/                                                                    _/
        !_/                   Math tools for GJK - EPA                         _/
        !_/                                                                    _/
        !_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
        !---------------------------------------------------------------
        !
        ! 3d�������
        PURE FUNCTION CROSS_PRODUCT_3D(A_, B_) RESULT(res_)
            IMPLICIT NONE
            REAL*8, INTENT(IN) :: A_(3), B_(3)
            REAL*8 :: res_(3) ! return value
            !  ->  ->  |  i   j   k |
            !  a X b = | a1  a2  a3 |
            !          | b1  b2  b3 |
            res_ = [A_(2) * B_(3) - A_(3) * B_(2), &
                    A_(3) * B_(1) - A_(1) * B_(3), &
                    A_(1) * B_(2) - A_(2) * B_(1)]
            RETURN
        END FUNCTION CROSS_PRODUCT_3D
        
        !---------------------------------------------------------------
        !
        ! �жϵ��Ƿ��ڵ������ڲ�
        FUNCTION isPointInSimplex(p_, simplex_) RESULT(res_)
            IMPLICIT NONE
            REAL*8, INTENT(IN) :: p_(3)
            REAL*8, INTENT(IN) :: simplex_(4,3)
            LOGICAL*1 :: res_ ! �������ʾԭ���Ƿ��ڵ������ڲ����߼�����
            LOGICAL*1 :: is ! 
            REAL*8 :: M(3)                     ! �����ε�����ƽ��
            REAL*8 :: MO(3)                     ! �����ε�����ƽ�� -> ԭ��
            REAL*8 :: AB(3), BC(3), dist(4), nml(4,3), vertexOnPlane(3,3)
            INTEGER*4 :: i, j ! 
            INTEGER*4, PARAMETER :: idFc(4, 3) = [1, 1, 1, 2, &
                                                  3, 2, 2 ,3, &
                                                  4, 4, 3, 4] ! ��������
            !---------------------------------------------------------------
            ! ��SITU 1���������ڵ����
            FORALL(i = 1:3) M(i) = SUM( simplex_(:,i) ) / 4.D0      

            ! ����������4���泯��ĵ�λ��ʸ
            ! face1 : [1 3 4],   face2 :[1 2 4]  face3 : [1 2 3]  face4 : [2 3 4]
            DO i = 1, 4, 1 
                AB = simplex_(idFc(i, 1), :) - simplex_(idFc(i, 2), :)
                BC = simplex_(idFc(i, 2), :) - simplex_(idFc(i, 3), :)
                nml(i,:) = UTZVEC( CROSS_PRODUCT_3D( AB, BC ) )
                IF ( DOT_PRODUCT( nml(i,:), simplex_(i,:) - M  ) < .0D0 ) nml(i,:) = - nml(i,:)
            END DO
            
            ! �����p_��������ľ���
            FORALL( i = 1:4) dist(i) = DOT_PRODUCT(simplex_(i,:) - p_, nml(i,:))
            
            ! ��SITU 2������������ϵ����
            DO i = 1, 4, 1
                IF ( DABS(dist(i) ) < 1.D-8 ) THEN
                    ! ��ȡ������Ƭ�ϵĶ�������
                    FORALL(j = 1:3) vertexOnPlane(j, :) = simplex_(idFc(i, j), :)
                    ! ���뺯�����ж�
                    is = IS_INSIDE_PF(vertexOnPlane, p_)
                    
                    IF ( is ) THEN
                        res_ = .TRUE.
                        RETURN
                    END IF
                      
                END IF
            END DO
            
            res_ = MERGE(.TRUE., .FALSE., ALL(dist > 0.D0))
            
            RETURN
        END FUNCTION isPointInSimplex
        
        
        !---------------------------------------------------------------
        !
        ! ����ƽ�������ϵ��������㷨���������ڶ���α��ϣ�
        PURE FUNCTION IS_INSIDE_PF(vertexOnPlane_, arbitraryPoint_) RESULT(res_)
            IMPLICIT NONE
            LOGICAL*1 :: res_ ! 
            REAL*8, INTENT(IN) :: vertexOnPlane_( :, :)
            REAL*8, INTENT(IN) :: arbitraryPoint_(3)
            !---------------------------------------------------------------
            REAL*8 :: temp
            INTEGER*4 :: i, NNODE
            REAL*8 :: crossProdResu( SIZE(vertexOnPlane_, 1) )    ! �洢��˽���������ţ�
            REAL*8 :: V( SIZE(vertexOnPlane_, 1) , 3)
            LOGICAL*1 :: zeroMask( SIZE(vertexOnPlane_, 1) )
            !---------------------------------------------------------------
            V = vertexOnPlane_
            zeroMask = .FALSE. 
            NNODE = SIZE(vertexOnPlane_, 1)
            !---------------------------------------------------------------
            ! ! �����  �ʵ�-�ǵ� X ������
            !����ͶӰ�� XOYƽ�棬���һ��ƽ������
            DO i = 1, NNODE, 1
                IF ( i == NNODE) THEN
                    crossProdResu(i) = (V(1,1) - V(i,1)) * (arbitraryPoint_(2) - V(i,2)) &
                                    - (V(1,2) - V(i,2)) * (arbitraryPoint_(1) - V(i,1))
                    EXIT
                END IF

                crossProdResu(i) = (V(i+1,1) - V(i,1)) * (arbitraryPoint_(2) - V(i,2)) &
                                   - (V(i+1,2) - V(i,2)) * (arbitraryPoint_(1) - V(i,1))
            END DO

            ! ����
            FORALL ( i = 1 : NNODE, DABS(crossProdResu(i)) < 1.0D-12 ) crossProdResu(i) = .0D0
            
            ! Ϊ�˱��ⱻ�жϵ������ʾ�ĵ�Ͷ��XOYƽ�����ֶ�㹲�߶�����crossProdResu=0�����
            ! ������������crossProdResu = 0 ������Ͷ��XOZƽ��������һ���ж�
            DO i = 1, NNODE, 1
                IF ( crossProdResu(i) > 1.0D-15 ) THEN ! ���ַ���
                    zeroMask(i) = .TRUE.
                END IF
            END DO
            IF ( ANY(zeroMask) == .FALSE. ) THEN !ȫΪ0
                !����ͶӰ�� XOZƽ�棬���һ��ƽ������
                DO i = 1, NNODE, 1
                    IF ( i == NNODE) THEN
                        crossProdResu(i) = (V(1,1) - V(i,1)) * (arbitraryPoint_(3) - V(i,3)) &
                                        - (V(1,3) - V(i,3)) * (arbitraryPoint_(1) - V(i,1))
                        EXIT
                    END IF

                    crossProdResu(i) = (V(i+1,1) - V(i,1)) * (arbitraryPoint_(3) - V(i,3)) &
                                       - (V(i+1,3) - V(i,3)) * (arbitraryPoint_(1) - V(i,1))
                END DO
            END IF
            
            
            ! ��������ÿһ��Ԫ�ض����һ��Ԫ����ˣ� �����ָ��ģ�˵�����ֲ�ͬ��
            DO i = 1, NNODE, 1
                temp = crossProdResu(1) * crossProdResu(i)
                ! �����ָ��ģ�˵�����ֲ�ͬ�ţ���ô�õ㲻��ƽ����
                IF ( temp < .0D0 ) THEN
                    res_ = .FALSE.
                    RETURN
                END IF
            END DO

            res_ = .TRUE.
            RETURN
        END FUNCTION IS_INSIDE_PF
        
        
        !---------------------------------------------------------------
        !
        ! ��ȡ��λ����
        PURE FUNCTION UTZVEC(vtr_) RESULT(res_)
            IMPLICIT NONE
            REAL*8, INTENT(IN) :: vtr_(:)
            REAL*8, ALLOCATABLE :: res_(:) ! return value
            REAL*8 :: MD
            ALLOCATE( res_( SIZE(vtr_) ) )
            MD = NORM2( vtr_ )
            res_ = MERGE(.0D0, vtr_ / MD, MD < 1.D-12)
            RETURN
        END FUNCTION UTZVEC
        
        !---------------------------------------------------------------
        !
        ! �㵽����ƽ��Ĵ�ֱ����(�����ţ�����ж�)
        FUNCTION DIST_PF_SIGN(arbitraryPoint_, defi3PoinOnPlane_) RESULT(res_)
            IMPLICIT NONE
            REAL*8 :: res_ ! return value
            REAL*8, INTENT(IN) :: arbitraryPoint_(3), defi3PoinOnPlane_(3,3)
            !---------------------------------------------------------------
            REAL*8 :: x_xyz(3) , p_xyz(3) 
            REAL*8 :: n_vtr(3) 
            x_xyz = arbitraryPoint_
            p_xyz = defi3PoinOnPlane_(1,:)
            n_vtr = UNINML(defi3PoinOnPlane_)
            
            ! �������ĵ㲻��ȷ��һ���棬��ô����
            IF ( ALL( DABS(n_vtr) < 1.D-12 ) ) THEN
                WRITE(UNIT=6, FMT="(A)") "ERROR - PURE FUNCTION DIST_PF_SIGN(arbitraryPoint_, defi3PoinOnPlane_) RESULT(res_)" ! 
                READ(UNIT=5, FMT=*)   
                STOP
            END IF
            
            res_ = DOT_PRODUCT(x_xyz - p_xyz, n_vtr)            
            RETURN                
        END FUNCTION DIST_PF_SIGN
        
        !---------------------------------------------------------------
        !
        ! ����ƽ��ĵ�λ��ʸ
        PURE FUNCTION UNINML(plane3Vertex_) RESULT(res_)
            IMPLICIT NONE
            REAL*8, INTENT(IN) :: plane3Vertex_(3,3)
            REAL*8 :: res_(3)
            !---------------------------------------------------------------
            REAL*8 :: N1(3), N2(3), cross(3)
            !---------------------------------------------------------------
            N1 = plane3Vertex_(2,:) - plane3Vertex_(1,:)
            N2 = plane3Vertex_(3,:) - plane3Vertex_(2,:)
            cross = CROSS_PRODUCT_3D(N1, N2)
            res_ = MERGE(cross / NORM2(cross), .0D0, ANY( DABS(cross) > 1.D-12 ))
            RETURN
        END FUNCTION UNINML
        
        !---------------------------------------------------------------
        !
        ! �жϸ�����һϵ�пռ���Ƿ�ȫ���غ� 
        PURE FUNCTION OVERLAP(a_) RESULT(res_)
            IMPLICIT NONE
            REAL*8, DIMENSION(:, :), INTENT(IN) :: a_
            LOGICAL :: res_
            INTEGER :: i, j, n_points
            REAL*8, PARAMETER :: TOLERANCE = 1.D-12

            n_points = SIZE(a_, 1)
            res_ = .TRUE.

            DO i = 1, n_points - 1
            DO j = i + 1, n_points
                IF (ANY(ABS(a_(i, :) - a_(j, :)) > TOLERANCE)) THEN
                res_ = .FALSE.
                RETURN
                END IF
            END DO
            END DO
            RETURN
        END FUNCTION OVERLAP
        
        !---------------------------------------------------------------
        !
        ! ����ֱ��ĳһ��Ĵ���ָ��õ�ĵ�λ����
        PURE FUNCTION VEC_PL(arbitraryPoint_, defi2PoinOnLine_) RESULT(res_)
            IMPLICIT NONE
            REAL*8 :: res_(3) ! return value
            REAL*8, INTENT(IN) :: arbitraryPoint_(3), defi2PoinOnLine_(2,3)
            !---------------------------------------------------------------
            REAL*8 :: A(3), B(3), AB(3), C(3), D(3), AC(3)
            REAL*8 :: vec(3)    ! ����
            !---------------------------------------------------------------
            A = defi2PoinOnLine_(1,:)
            B = defi2PoinOnLine_(2,:)
            C = arbitraryPoint_
            AB = B - A
            AC = C - A
            D = A + DOT_PRODUCT(AC, AB) / NORM2(AB) * UTZVEC(AB)
            !---------------------------------------------------------------
            res_ = UTZVEC(D - C)
            RETURN
        END FUNCTION VEC_PL
        

        !---------------------------------------------------------------
        !
        ! ��ȡ�ڿռ��ϵõ���ֱ�߼���̾��룬��Ӧ�ֱ�������ֱ���ϵĵ㣨2�����㣩
        FUNCTION FOOT_LL(lineDefiEP_1_, lineDefiEP_2_) RESULT(res_)
            IMPLICIT NONE
            REAL*8, INTENT(IN) :: lineDefiEP_1_(2,3), lineDefiEP_2_(2,3)
            REAL*8 :: res_(2,3)
            !---------------------------------------------------------------
            REAL*8 :: P1(3) , Q1(3) , P2(3) , Q2(3) 
            REAL*8 :: d1_vec(3) , d2_vec(3) , r_vec(3) 
            REAL*8 :: a , b , c , d , e , f 
            REAL*8 :: s , t 
            REAL*8 :: L1_s(3) , L2_t(3) 
            !---------------------------------------------------------------
            P1 = lineDefiEP_1_(1,:)
            Q1 = lineDefiEP_1_(2,:)
            P2 = lineDefiEP_2_(1,:)
            Q2 = lineDefiEP_2_(2,:)

            d1_vec = Q1 - P1
            d2_vec = Q2 - P2
            r_vec = P1 - P2

            a = DOT_PRODUCT(d1_vec, d1_vec)
            b = DOT_PRODUCT(d1_vec, d2_vec)
            c = DOT_PRODUCT(d1_vec, r_vec)
            e = DOT_PRODUCT(d2_vec, d2_vec)
            f = DOT_PRODUCT(d2_vec, r_vec)

            d = a * e - b**2

            IF ( DABS(d) < 1.0D-12 ) THEN !��ֱ��ƽ�У���ôȡ��һ���ߵ��е���Ϊ����
                res_(1,:) = (P1 + Q1) / 2.D0 
                res_(2,:) = FOOT_PL(res_(1,:), lineDefiEP_2_)
            ELSE !��ֱ�߲�ƽ��
                s = (b * f - c * e) / d
                t = (a * f - b * c) / d
                L1_s = P1 + s * (Q1 - P1)
                L2_t = P2 + t * (Q2 - P2)

                res_(1,:) = L1_s
                res_(2,:) = L2_t
            END IF
            RETURN
        END FUNCTION FOOT_LL
        
        !---------------------------------------------------------------
        !
        ! �㵽����ֱ�ߵĴ���
        PURE FUNCTION FOOT_PL(arbitraryPoint_, defi2PoinOnLine_) RESULT(res_)
            IMPLICIT NONE
            REAL*8 :: res_(3) ! return value
            REAL*8, INTENT(IN) :: arbitraryPoint_(3), defi2PoinOnLine_(2,3)
            !---------------------------------------------------------------
            REAL*8 :: P(3)   ! ���жϵĵ�
            REAL*8 :: V(2,3)  ! ���жϵ�ֱ���ϵ���������
            !---------------------------------------------------------------
            P = arbitraryPoint_
            V = defi2PoinOnLine_
            !---------------------------------------------------------------
            res_ = V(1,:) + DOT_PRODUCT( P - V(1,:), UTZVEC(V(2,:) - V(1,:)) ) * UTZVEC(V(2,:) - V(1,:))
            RETURN
        END FUNCTION FOOT_PL
        
        !---------------------------------------------------------------
        !
        ! ���ռ���ɢ�ҵĵ㰴ʱ�����У����������㶼��ͬһ�ռ�ƽ���ڣ����������ɿ�
        !  IN - points(:,3)
        !  OUT - ordered_points(:,3)
        !---------------------------------------------------------------
        PURE FUNCTION SORT_CLOCK(points_) RESULT(ordered_points_)
            REAL*8, DIMENSION(:, :), INTENT(IN) :: points_
            REAL*8, DIMENSION(SIZE(points_, 1), 3) :: ordered_points_
            REAL*8, DIMENSION(3) :: centroid, normal, v1, v2
            INTEGER*4 :: i, j, index, num_points
            REAL*8 :: angle, min_angle

            ! ���жϼ������Ƿ��غϣ��غ���ֱ���������㷨����������
            IF ( OVERLAP(points_) ) RETURN
            
            num_points = SIZE(points_, 1)

            ! ��������
            centroid = SUM(points_, DIM=1) / num_points

            ! ����ƽ�淨����
            v1 = points_(2, :) - points_(1, :)
            v2 = points_(3, :) - points_(1, :)
            normal = CROSS_PRODUCT_3D(v1, v2)

            ! �Ե�һ������Ϊ��׼�����Ƕ����������
            ordered_points_(1, :) = points_(1, :)
            DO i = 2, num_points
                min_angle = HUGE(angle)
                index = -1
                DO j = 1, num_points
                    ! ����������ĵ�
                    IF (IS_POINT_IN_ORDERED_POINTS(points_(j, :), ordered_points_)) CYCLE
                    v1 = points_(j, :) - centroid
                    v2 = ordered_points_(i - 1, :) - centroid
                    ! CHOOSE : ��ʱ������
                    angle = ATAN2(DOT_PRODUCT(normal, CROSS_PRODUCT_3D(v2, v1)), DOT_PRODUCT(v1, v2))
                    ! CHOOSE : ˳ʱ������
                    !angle = atan2(dot_product(normal, cross_product_for_fn_order_points(v1, v2)), dot_product(v1, v2))
                    angle = MODULO(angle + 2.0 * ACOS(-1.0), 2.0 * ACOS(-1.0))

                    IF (angle < min_angle) THEN
                        min_angle = angle
                        index = j
                    END IF
                END DO
                ordered_points_(i, :) = points_(index, :)
            END DO
            RETURN
            CONTAINS
                !---------------------------------------------------------------
                ! �����Ƿ��Ѿ��������ĵ㼯��
                PURE FUNCTION IS_POINT_IN_ORDERED_POINTS(point, ordered_points) RESULT(is_in)
                    REAL*8, DIMENSION(3), INTENT(IN) :: point
                    REAL*8, DIMENSION(:, :), INTENT(IN) :: ordered_points
                    LOGICAL*1 :: is_in
                    INTEGER*4 :: i
                    is_in = .FALSE.
                    DO i = 1, SIZE(ordered_points, 1)
                        IF (ALL(point == ordered_points(i, :))) THEN
                            is_in = .TRUE.
                            EXIT
                        END IF
                    END DO
                    RETURN
                END FUNCTION IS_POINT_IN_ORDERED_POINTS
                
        END FUNCTION SORT_CLOCK
            
            
        PURE FUNCTION GET_RANDOM_UNIT_VECTOR(index_) RESULT(res_)
            IMPLICIT NONE
            INTEGER*4, INTENT(IN) :: index_ ! 
            REAL*8 :: res_(3)
            REAL*8, PARAMETER :: dataLib(3,100) = &
            RESHAPE(&
             [ 0.000001109357820885D0,  0.072093544214837393D0,  0.997397874913172555D0, &
               0.266483497218669374D0, -0.727347325988231153D0,  0.632417910157418883D0, &
               0.079214616132658941D0, -0.782543920607548071D0, -0.617535470164364719D0, &
              -0.993301267605208316D0,  0.106810772229378015D0,  0.044091390425458579D0, &
               0.082261341377368513D0,  0.991595302008176138D0, -0.099859044408155587D0, &
              -0.787452696781838490D0,  0.616178410256023601D0,  0.015569748404171571D0, &
              -0.247966562512464128D0,  0.750010049461640738D0, -0.613186357955148420D0, &
              -0.715817591888975313D0,  0.423804523888427931D0,  0.554972882827594716D0, &
               0.499764308041154848D0,  0.237809719054367125D0, -0.832875845448425078D0, &
               0.360748686617363812D0,  0.307777557994801998D0,  0.880416583157429655D0, &
               0.713138609686784886D0, -0.678418744074228530D0,  0.176582363396647901D0, &
               0.881992030996567422D0,  0.026379550968972942D0, -0.470525426039045791D0, &
              -0.267765386517834436D0,  0.464539693453386748D0, -0.844099858422679872D0, &
               0.513202226307113540D0,  0.794177664474205347D0,  0.325430963744568147D0, &
               0.266257765457365569D0,  0.689919118649417573D0,  0.673140707471819200D0, &
              -0.533214734590422568D0,  0.393416539739102400D0,  0.748936227642498564D0, &
              -0.623072641479377243D0, -0.654446770357797636D0,  0.428345547669355065D0, &
               0.584825748689458469D0,  0.437231667603634577D0,  0.683232985528625658D0, &
              -0.556342156780530561D0, -0.693940941632379182D0, -0.457087928209829908D0, &
               0.797251122953163582D0, -0.186816815361580540D0, -0.574012303394340728D0, &
               0.652717880922520921D0,  0.670487884243165855D0, -0.352711447230079078D0, &
              -0.119569576931363289D0, -0.933186657472575787D0, -0.338918542702544345D0, &
               0.662896092871913201D0, -0.734670864402726664D0,  0.144317327625279795D0, &
              -0.453865743569666802D0,  0.555714019359183631D0,  0.696554244478931106D0, &
               0.654083844194692787D0, -0.209153829113278511D0,  0.726931221320659904D0, &
               0.590510679076412859D0,  0.337909209878702432D0,  0.732880961531860775D0, &
               0.968625410428645917D0, -0.064469344047131227D0, -0.240017745073296679D0, &
              -0.836672384182689188D0, -0.337478629755403936D0,  0.431378599381644634D0, &
               0.415710848734430150D0,  0.722574771381445879D0, -0.552331594250728086D0, &
              -0.333326475889782536D0,  0.815058361243497620D0, -0.473891684077661635D0, &
              -0.652533192903382075D0, -0.591467557663984178D0,  0.473673474442383280D0, &
               0.394659527294562162D0, -0.550384256978558417D0, -0.735745218935055623D0, &
              -0.636304506189762753D0,  0.473703705794754570D0,  0.608868930492367122D0, &
              -0.719230459123433086D0, -0.158162890699728137D0,  0.676529413015133918D0, &
               0.629759138526492901D0, -0.491788561913722666D0,  0.601288148738358452D0, &
               0.584411917965700356D0, -0.367877772023600003D0,  0.723276333769192092D0, &
               0.870106618562407896D0, -0.204182999880998167D0,  0.448579730809907151D0, &
               0.529356795812083503D0, -0.718211329438827373D0,  0.451612520855297239D0, &
               0.733690094242977708D0, -0.622391387307088984D0, -0.272631264926984196D0, &
              -0.605777076602946218D0, -0.315061533953294726D0,  0.730595896022818714D0, &
              -0.761009425976650333D0, -0.636619547995314727D0,  0.124820690131605891D0, &
              -0.646761961270369112D0, -0.761942845893679443D0, -0.033794452875378959D0, &
               0.365154502536077674D0,  0.505749055061637143D0, -0.781588179658502025D0, &
               0.574247267419540908D0,  0.634851826576257938D0,  0.516917047652695638D0, &
               0.346341716472641781D0, -0.550932683186623917D0, -0.759289532410115098D0, &
              -0.794685184986554050D0, -0.055389826957407198D0,  0.604489391000797349D0, &
              -0.416259521322270454D0, -0.054995592820233065D0, -0.907581123469910711D0, &
               0.794777927582307919D0,  0.342095783921817331D0, -0.501296838660377997D0, &
              -0.338337965454608924D0, -0.286035970801144568D0, -0.896499216140138389D0, &
              -0.726532004741409887D0, -0.049688151104356579D0, -0.685333738937649595D0, &
              -0.603734615736470803D0, -0.585014438414317439D0,  0.541537275363678683D0, &
              -0.676560375498003186D0, -0.722348934167962309D0,  0.143101626868494480D0, &
               0.586582880385051575D0,  0.072766280975167824D0, -0.806613657702508258D0, &
              -0.755532705527683479D0, -0.071266043707085253D0, -0.651223066155029895D0, &
              -0.920701606636518566D0,  0.311540070620156373D0,  0.235056027225258340D0, &
               0.541712171882508864D0, -0.838526306892959261D0,  0.058494063654270075D0, &
              -0.408115455093796653D0, -0.092597310866135374D0, -0.908222171791651101D0, &
              -0.258240219479359101D0, -0.908622337155473581D0,  0.328203347736395479D0, &
              -0.061612129227968819D0, -0.446992987857170232D0,  0.892413141061087156D0, &
               0.788042672316281223D0, -0.496244917147545261D0,  0.364320914598434853D0, &
              -0.248619129130190686D0,  0.619445212796131295D0, -0.744631557869058658D0, &
               0.727207891810358387D0, -0.392604991169558049D0, -0.563054174123134521D0, &
              -0.730052156895783066D0,  0.157234865285751285D0,  0.665057174497340808D0, &
               0.600414670664006778D0,  0.750265884008508910D0,  0.276773059643389052D0, &
              -0.083928500830154310D0,  0.690568080639724524D0,  0.718381328230327632D0, &
               0.694831042024156353D0,  0.584804220606428005D0, -0.418585530806468986D0, &
              -0.111848450943919986D0, -0.781531383436509852D0, -0.613757786692161189D0, &
              -0.279182094755242194D0, -0.930461735000781665D0, -0.237272665234346397D0, &
              -0.689964963785805074D0, -0.305025070889099192D0, -0.656435872631251471D0, &
               0.633382581384791088D0,  0.583236672149216373D0,  0.508587740570539015D0, &
               0.466924244038473768D0, -0.606103736912413371D0,  0.643909939688702027D0, &
              -0.137658227056735444D0, -0.193627586092586290D0, -0.971369430457616478D0, &
               0.393853240338342958D0,  0.768953844741995574D0,  0.503576816117948800D0, &
              -0.132535470218959284D0,  0.729368436809752718D0, -0.671160213748950629D0, &
               0.159029880166712406D0,  0.267247506574191773D0,  0.950414787050390064D0, &
               0.585440601303706010D0, -0.650059126571057910D0,  0.484440331007677694D0, &
               0.086766095195569742D0, -0.926700911609081412D0,  0.365646092755564367D0, &
               0.404761320436991479D0, -0.409969869053845359D0, -0.817369549191842681D0, &
              -0.630382450683336315D0,  0.770188809015893039D0, -0.097093585458315548D0, &
              -0.042053492941287379D0, -0.611271645428856480D0, -0.790302776931813389D0, &
               0.929725661108754209D0,  0.077330619173836948D0, -0.360041900914436386D0, &
              -0.889604251783720934D0, -0.344981229410519730D0, -0.299319606044663511D0, &
               0.129702915764274479D0, -0.696106796017660678D0, -0.706124976318124986D0, &
              -0.796994723739967381D0, -0.420325416758673909D0, -0.433734889485847597D0, &
              -0.643021987392653815D0, -0.525087908251825164D0,  0.557499248732520325D0, &
               0.223259530927500754D0, -0.439307839166757808D0,  0.870151598456651798D0, &
               0.639217882809690274D0,  0.671377686488942249D0,  0.375036665382270096D0, &
               0.228323372420344811D0, -0.748223967023273318D0, -0.622920005119883879D0, &
              -0.632452534964462632D0,  0.397443937197173747D0, -0.664862472848508856D0, &
              -0.575267651846246730D0,  0.586755089131675400D0,  0.569899635126559057D0, &
               0.934572561750450670D0,  0.355419405776895792D0,  0.015848432742659273D0, &
              -0.122211293462219608D0,  0.261591882966958789D0,  0.957410093176425669D0, &
               0.418206651287156450D0, -0.714638510825073237D0,  0.560709368269252773D0, &
              -0.455037020713617735D0,  0.389115382040291002D0,  0.800956009553404180D0, &
               0.576937065595787169D0, -0.543479726634975457D0,  0.609732243758270287D0, &
              -0.094516770591717383D0,  0.753943490941892613D0,  0.650104447410771891D0, &
               0.489068888565033721D0, -0.424755340422356520D0,  0.761836283607213560D0, &
               0.986861350764715373D0,  0.139794765568494128D0,  0.081006776793618909D0, &
              -0.902962972513389861D0, -0.262938852206923646D0,  0.339883848203895222D0, &
              -0.712980642840275625D0,  0.087812143183863101D0,  0.695663446247195227D0], &
            [3,100])
                
            res_ = dataLib(:, index_)
            
            RETURN
        END FUNCTION GET_RANDOM_UNIT_VECTOR
            
            
    END MODULE GCLIB_GJKEPA


