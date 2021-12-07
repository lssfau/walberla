target_sources( pe
    PRIVATE
    RigidBodyRemovalNotification.h
    RigidBodyCopyNotification.h
    RigidBodyRemoteMigrationNotification.h
    RigidBodyMigrationNotification.h
    RigidBodyVelocityCorrectionNotification.h
    RigidBodyNewShadowCopyNotification.h
    RigidBodyDeletionNotification.h
    ParseMessage.h
    PackNotification.h
    Instantiate.h
    RigidBodyVelocityUpdateNotification.h
    Marshalling.h
    RigidBodyUpdateNotification.h
    Marshalling.cpp
    RigidBodyForceNotification.h
    DynamicMarshalling.h
    RigidBodyRemovalInformationNotification.h
    NotificationType.h     
    )
add_subdirectory( rigidbody )